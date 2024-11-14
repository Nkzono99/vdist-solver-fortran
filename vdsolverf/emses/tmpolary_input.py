from pathlib import Path
from typing import Any, List

import emout
import f90nml
import f90nml.namelist
import numpy as np

from .geotype import (
    create_cylinder_boundary,
    create_rectangular_boundary,
    create_sphere_boundary,
)

TMP_INP_KEYS = {
    "esorem": ["emflag"],
    "plasma": ["wp", "wc", "phixy", "phiz"],
    "tmgrid": ["dt", "nx", "ny", "nz"],
    "system": ["nspec", "npbnd"],
    "intp": ["qm", "path", "peth", "vdri", "vdthz", "vdthxy", "spa", "spe", "speth"],
    "ptcond": [
        "zssurf",
        "xlrechole",
        "xurechole",
        "ylrechole",
        "yurechole",
        "zlrechole",
        "zurechole",
        "boundary_type",
        "boundary_types",
        "cylinder_origin",
        "cylinder_radius",
        "cylinder_height",
        "rcurv",
        "rectangle_shape",
        "sphere_origin",
        "sphere_radius",
        "circle_origin",
        "circle_radius",
        "cuboid_shape",
        "disk_origin",
        "disk_height",
        "disk_radius",
        "disk_inner_raidus",
        "plane_with_circle_hole_zlower",
        "plane_with_circle_hole_height",
        "plane_with_circle_hole_radius",
    ],
    "emissn": [
        "nflag_emit",
        "nepl",
        "curf",
        "nemd",
        "curfs",
        "xmine",
        "xmaxe",
        "ymine",
        "ymaxe",
        "zmine",
        "zmaxe",
        "thetaz",
        "thetaxy",
    ],
}


class TempolaryInput(object):
    def __init__(self, data: emout.Emout):
        self.__data = data
        self.__tmppath: Path = data.directory / f"plasma-vdsolverf.inp"

    def __enter__(self) -> "TempolaryInput":
        data = self.__data
        from collections import defaultdict

        dic = defaultdict(lambda: dict())

        for group in TMP_INP_KEYS.keys():
            for key in TMP_INP_KEYS[group]:
                if key not in data.inp:
                    continue
                dic[group][key] = getattr(data.inp, key)

        inp = f90nml.Namelist(dic)

        # start_indexを追加する.
        for group in TMP_INP_KEYS.keys():
            start_indexes = data.inp.nml[group].start_index
            for key in TMP_INP_KEYS[group]:
                if key not in start_indexes:
                    continue
                inp[group].start_index[key] = start_indexes[key]

        self.convert_from_geotype(inp)

        inp.write(str(self.__tmppath.resolve()), force=True)

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def convert_from_geotype(self, nml: f90nml.Namelist):
        data = self.__data

        if "geotype" not in data.inp:
            return

        if "boundary_type" in data.inp and data.inp.boundary_type != "complex":
            nml["ptcond"]["boundary_types"] = [data.inp.boundary_type]
            nml["ptcond"].start_index["boundary_types"] = [1]

        if "boundary_type" not in data.inp:
            nml["ptcond"]["boundary_type"] = "complex"
            nml["ptcond"]["boundary_types"] = []
            nml["ptcond"].start_index["boundary_types"] = [1]

        if "npc" not in data.inp:
            return

        for ipc in range(data.inp.npc):
            if data.inp.geotype[ipc] in (0, 1):
                create_rectangular_boundary(nml, data, ipc)
            elif data.inp.geotype[ipc] == 2:
                create_cylinder_boundary(nml, data, ipc)
            elif data.inp.geotype[ipc] == 3:
                create_sphere_boundary(nml, data, ipc)
            else:
                raise NotImplementedError()

    @property
    def tmppath(self):
        return self.__tmppath
