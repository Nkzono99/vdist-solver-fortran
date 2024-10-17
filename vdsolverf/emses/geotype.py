import emout
import f90nml
import f90nml.namelist
import numpy as np
from typing import List


def create_rectangular_boundary(
    nml: f90nml.Namelist, data: emout.Emout, ipc: int
) -> None:
    xlpc = fetch_from_inp(data, "ptcond", "xlpc", default=0.0)[ipc]
    ylpc = fetch_from_inp(data, "ptcond", "ylpc", default=0.0)[ipc]
    zlpc = fetch_from_inp(data, "ptcond", "zlpc", default=0.0)[ipc]
    xupc = fetch_from_inp(data, "ptcond", "xupc", default=0.0)[ipc]
    yupc = fetch_from_inp(data, "ptcond", "yupc", default=0.0)[ipc]
    zupc = fetch_from_inp(data, "ptcond", "zupc", default=0.0)[ipc]

    append_cuboid_to_namelist(nml, [xlpc, xupc, ylpc, yupc, zlpc, zupc])


def create_cylinder_boundary(nml: f90nml.Namelist, data: emout.Emout, ipc: int) -> None:
    bdyalign = fetch_from_inp(data, "ptcond", "bdyalign", default=0)
    bdyedge = fetch_from_inp(data, "ptcond", "bdyedge", default=0)
    bdyradius = fetch_from_inp(data, "ptcond", "bdyradius", default=0)
    bdycoord = fetch_from_inp(data, "ptcond", "bdycoord", default=0)

    height = bdyedge[ipc, 1] - bdyedge[ipc, 0]
    radius = bdyradius[ipc]

    if bdyalign[ipc] == 1:
        axis = "x"
        origin = np.array([bdyedge[ipc, 0], bdycoord[ipc, 0], bdycoord[ipc, 1]])
        upper_origin = origin + np.array([height, 0, 0])
    elif bdyalign[ipc] == 2:
        axis = "y"
        origin = np.array([bdycoord[ipc, 1], bdyedge[ipc, 0], bdycoord[ipc, 0]])
        upper_origin = origin + np.array([0, height, 0])
    else:
        axis = "z"
        origin = np.array([bdycoord[ipc, 0], bdycoord[ipc, 1], bdyedge[ipc, 0]])
        upper_origin = origin + np.array([0, 0, height])

    origin = origin.tolist()
    upper_origin = upper_origin.tolist()

    append_circle_to_namelist(nml, axis, origin, radius)
    append_circle_to_namelist(nml, axis, upper_origin, radius)
    append_cylinder_to_namelist(nml, axis, origin, height, radius)


def create_sphere_boundary(nml: f90nml.Namelist, data: emout.Emout, ipc: int) -> None:
    bdycoord = fetch_from_inp(data, "ptcond", "bdycoord", default=0.0)
    bdyradius = fetch_from_inp(data, "ptcond", "bdyradius", default=0.0)

    origin = bdycoord[ipc, :].tolist()
    radius = bdyradius[ipc]

    append_sphere_to_namelist(nml, origin, radius)


def append_cuboid_to_namelist(nml: f90nml.namelist, shape):
    index = len(nml["ptcond"]["boundary_types"])

    append_to_namelist(
        nml=nml,
        group="ptcond",
        name="boundary_types",
        value="cuboid",
        index=index,
        ndim=1,
    )

    append_to_namelist(
        nml=nml, group="ptcond", name="cuboid_shape", value=shape, index=index, ndim=2
    )


def append_circle_to_namelist(nml, axis, origin, radius):
    index = len(nml["ptcond"]["boundary_types"])

    append_to_namelist(
        nml=nml,
        group="ptcond",
        name="boundary_types",
        value=f"circle{axis}",
        index=index,
        ndim=1,
    )
    append_to_namelist(
        nml=nml, group="ptcond", name="circle_origin", value=origin, index=index, ndim=2
    )
    append_to_namelist(
        nml=nml, group="ptcond", name="circle_radius", value=radius, index=index, ndim=1
    )


def append_cylinder_to_namelist(
    nml, axis: str, origin: List[float], radius: float, height: float
):
    index = len(nml["ptcond"]["boundary_types"])

    append_to_namelist(
        nml=nml,
        group="ptcond",
        name="boundary_types",
        value=f"cylinder{axis}",
        index=index,
        ndim=1,
    )
    append_to_namelist(
        nml=nml,
        group="ptcond",
        name="cylinder_origin",
        value=origin,
        index=index,
        ndim=2,
    )
    append_to_namelist(
        nml=nml,
        group="ptcond",
        name="cylinder_radius",
        value=radius,
        index=index,
        ndim=1,
    )
    append_to_namelist(
        nml=nml,
        group="ptcond",
        name="cylinder_height",
        value=height,
        index=index,
        ndim=1,
    )


def append_sphere_to_namelist(nml, origin, radius):
    index = len(nml["ptcond"]["boundary_types"])

    append_to_namelist(
        nml=nml,
        group="ptcond",
        name="boundary_types",
        value=f"sphere",
        index=index,
        ndim=1,
    )
    append_to_namelist(
        nml=nml, group="ptcond", name="sphere_origin", value=origin, index=index, ndim=2
    )
    append_to_namelist(
        nml=nml, group="ptcond", name="sphere_radius", value=radius, index=index, ndim=1
    )


def append_to_namelist(
    nml: f90nml.Namelist, group: str, name: str, value, index: int, ndim: int
):
    if name not in nml[group]:
        nml[group][name] = [value]
        nml[group].start_index[name] = [1] * (ndim - 1) + [index + 1]
        return

    start_index = nml[group].start_index[name][-1]
    for _ in range(len(nml[group][name]) + start_index - 1, index):
        none_value = None
        for _ in range(ndim - 1):
            none_value = [none_value] * 3
        nml[group][name].append(none_value)

    nml[group][name].append(value)


def fetch_from_inp(data: emout.Emout, group: str, name: str, default=0) -> np.ndarray:
    start_index = np.array(data.inp.nml[group].start_index[name][::-1]) - 1

    values = getattr(data.inp, name)
    values = np.array(values, dtype=float)
    values = np.nan_to_num(values, nan=default)

    original_shape = start_index + np.array(values.shape)

    original_values = np.full(original_shape, default)
    original_values[tuple(slice(i, None) for i in start_index)] = values

    return original_values
