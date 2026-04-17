from pathlib import Path

import emout
import f90nml

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
        "boundary_conductor_id",
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
        "disk_inner_radius",
        "conductivity",
        "plane_with_circle_hole_zlower",
        "plane_with_circle_hole_height",
        "plane_with_circle_hole_radius",
        "plane_with_circle_origin",
        "plane_with_circle_radius",
        "max_bounce_count",
        "boundary_mirror_reflection_rate",
        "boundary_mirror_reflect_alpha",
        "boundary_reversal_reflection_rate",
        "boundary_reversal_reflect_alpha",
        "boundary_mirror_reflect_energy_loss_frac",
        "boundary_reversal_reflect_energy_loss_frac",
        "enable_secondary_electron_emission",
        "boundary_se_model_type",
        "boundary_se_const_yield",
        "boundary_se_yield_max",
        "boundary_se_energy_max",
        "boundary_se_species_id",
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
        inp = f90nml.Namelist()

        for group, keys in TMP_INP_KEYS.items():
            group_namelist = self._create_group_namelist(group, keys)
            if group_namelist:
                inp[group] = group_namelist

        self.convert_from_geotype(inp)

        inp.write(str(self.__tmppath.resolve()), force=True)

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.__tmppath.exists():
            self.__tmppath.unlink()

    def _create_group_namelist(self, group: str, keys):
        data = self.__data

        if group not in data.inp.nml:
            return None

        source_group = data.inp.nml[group]
        group_namelist = f90nml.Namelist()

        for key in keys:
            if key not in data.inp:
                continue

            group_namelist[key] = getattr(data.inp, key)
            if key in source_group.start_index:
                group_namelist.start_index[key] = source_group.start_index[key]

        if not group_namelist:
            return None

        return group_namelist

    def _ensure_complex_boundary_group(self, nml: f90nml.Namelist) -> f90nml.Namelist:
        if "ptcond" not in nml:
            nml["ptcond"] = f90nml.Namelist()

        ptcond = nml["ptcond"]
        original_boundary_type = ptcond.get("boundary_type", "complex")

        existing_boundary_types = ptcond.get("boundary_types", [])
        if isinstance(existing_boundary_types, str):
            boundary_types = [existing_boundary_types]
        else:
            boundary_types = list(existing_boundary_types)

        if original_boundary_type != "complex" and original_boundary_type not in boundary_types:
            boundary_types.insert(0, original_boundary_type)

        ptcond["boundary_type"] = "complex"
        ptcond["boundary_types"] = boundary_types
        ptcond.start_index["boundary_types"] = [1]

        return ptcond

    def convert_from_geotype(self, nml: f90nml.Namelist):
        data = self.__data

        if "geotype" not in data.inp:
            return

        self._ensure_complex_boundary_group(nml)

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
