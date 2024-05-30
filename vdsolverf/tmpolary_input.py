from pathlib import Path

import emout
import f90nml

TMP_INP_KEYS = {
    "esorem": ["emflag"],
    "plasma": ["wp", "wc", "phixy", "phiz"],
    "tmgrid": ["dt", "nx", "ny", "nz"],
    "system": ["nspec", "npbnd"],
    "intp": ["qm", "path", "peth", "vdri", "vdthz", "vdthxy"],
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
        "circle_raidus",
        "cuboid_shape",
        "disk_origin",
        "disk_height",
        "disk_radius",
        "disk_inner_raidus",
    ],
    "emissn": [
        "nflag_emit",
        "nepl",
        "curf",
        "nemd",
        "curfs",
        "xmine",
        "xmaxe" "ymine",
        "ymaxe" "zmine",
        "zmaxe",
        "thetaz",
        "thetaxy",
    ],
}


class TempolaryInput(object):
    def __init__(self, data: emout.Emout):
        self.__data = data
        self.__tmppath: Path = data.directory / f"plasma.inp.{hash(str(data.inp))}"

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
        inp.write(str(self.__tmppath.resolve()))

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.__tmppath.unlink()

    @property
    def tmppath(self):
        return self.__tmppath
