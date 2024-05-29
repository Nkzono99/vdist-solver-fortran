from ctypes import *
from os import PathLike
from pathlib import Path
from typing import List, Tuple, Literal

import emout
import numpy as np
from vdsolver.core import Particle

from .tmpolary_input import TempolaryInput


VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX = (
    Path().home() / ".local/lib" / "libvdist-solver-fortran.so"
)


def get_probabirities(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: List[Particle],
    dt_multiplier: float,
    max_step: int,
    max_probabirity_types: int = 100,
    os: Literal["linux"] = "linux",
    library_path: PathLike = None,
):
    if os == "linux":
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX
        return get_probabirities_linux(
            directory=directory,
            ispec=ispec,
            istep=istep,
            particles=particles,
            dt_multiplier=dt_multiplier,
            max_step=max_step,
            max_probabirity_types=max_probabirity_types,
            library_path=library_path,
        )


def get_probabirities_linux(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: List[Particle],
    dt_multiplier: float,
    max_step: int,
    max_probabirity_types: int = 100,
    library_path: PathLike = VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX,
) -> Tuple[np.ndarray, List[Particle]]:
    dll = CDLL(library_path)

    dll.get_probabirities.argtypes = [
        c_char_p,  # inppath
        c_int,  # length
        c_int,  # lx
        c_int,  # ly
        c_int,  # lz
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),  # ebvalues
        c_int,  # ispec
        c_int,  # npcls
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # positions
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # velocities
        c_double,  # dt_multiplier
        c_int,  # max_step
        c_int,  # max_probabirity_types
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),  # return_probabirities
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # return_positions
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # return_velocities
    ]
    dll.get_probabirities.restype = None

    data = emout.Emout(directory)

    ebvalues = create_relocated_ebvalues(data, istep)

    return_probabirities = np.empty(len(particles), dtype=np.float64)
    return_positions = np.empty((len(particles), 3), dtype=np.float64)
    return_velocities = np.empty((len(particles), 3), dtype=np.float64)

    positions = np.array([particle.pos for particle in particles], dtype=np.float64)
    velocities = np.array([particle.vel for particle in particles], dtype=np.float64)

    with TempolaryInput(data) as tmpinp:
        inppath = tmpinp.tmppath
        inppath_str = str(inppath.resolve())

        _inppath = create_string_buffer(inppath_str.encode())
        _length = c_int(len(inppath_str))
        _nx = c_int(data.inp.nx)
        _ny = c_int(data.inp.ny)
        _nz = c_int(data.inp.nz)
        _ispec = c_int(ispec + 1)
        _nparticles = c_int(len(particles))
        _dt_multiplier = c_double(dt_multiplier)
        _max_step = c_int(max_step)
        _max_probabirity_types = c_int(max_probabirity_types)

        dll.get_probabirities(
            _inppath,
            _length,
            _nx,
            _ny,
            _nz,
            ebvalues,
            _ispec,
            _nparticles,
            positions,
            velocities,
            _dt_multiplier,
            _max_step,
            _max_probabirity_types,
            return_probabirities,
            return_positions,
            return_velocities,
        )

    handle = dll._handle
    cdll.LoadLibrary("libdl.so").dlclose(handle)

    return_particles = [
        Particle(pos, vel) for pos, vel in zip(return_positions, return_velocities)
    ]

    return_probabirities[return_probabirities == -1] = np.nan

    return return_probabirities, return_particles


def create_relocated_ebvalues(data, istep):
    ebvalues = np.zeros(
        (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1, 6), dtype=np.float64
    )

    # EX
    ielem = 0
    ex = data.ex[istep]
    ebvalues[:, :, 1:-1, ielem] = 0.5 * (ex[:, :, :-2] + ex[:, :, 1:-1])
    if data.inp.mtd_vbnd[0] in [0, 2]:
        ebvalues[:, :, 0, ielem] = 0
        ebvalues[:, :, -1, ielem] = 0
    else:
        ebvalues[:, :, 0, ielem] = ex[:, :, 0]
        ebvalues[:, :, -1, ielem] = ex[:, :, -1]

    # EY
    ielem = 1
    ey = data.ey[istep]
    ebvalues[:, 1:-1, :, ielem] = 0.5 * (ey[:, :-2, :] + ey[:, 1:-1, :])
    if data.inp.mtd_vbnd[0] in [0, 2]:
        ebvalues[:, 0, :, ielem] = 0
        ebvalues[:, -1, :, ielem] = 0
    else:
        ebvalues[:, 0, :, ielem] = ey[:, 0, :]
        ebvalues[:, -1, :, ielem] = ey[:, -1, :]

    # EZ
    ielem = 2
    ez = data.ez[istep]
    ebvalues[1:-1, :, :, ielem] = 0.5 * (ez[:-2, :, :] + ez[1:-1, :, :])
    if data.inp.mtd_vbnd[0] in [0, 2]:
        ebvalues[0, :, :, ielem] = 0
        ebvalues[-1, :, :, ielem] = 0
    else:
        ebvalues[0, :, :, ielem] = ez[0, :, :]
        ebvalues[-1, :, :, ielem] = ez[-1, :, :]

    # BX
    ielem = 3
    bx = data.bx[istep]
    # yz平面に1グリッド覆うように拡張する
    bxe = np.empty(np.array(bx.shape) + np.array([1, 1, 0]))
    bxe[1:-1, 1:-1, :] = bx[:-1, :-1, :]
    if data.inp.mtd_vbnd[1] == 0:
        bxe[1:-1, 0, :] = bxe[1:-1, -2, :]
        bxe[1:-1, -1, :] = bxe[1:-1, 1, :]
    elif data.inp.mtd_vbnd[1] == 1:
        bxe[1:-1, 0, :] = -bxe[1:-1, 1, :]
        bxe[1:-1, -1, :] = -bxe[1:-1, -2, :]
    else:
        bxe[1:-1, 0, :] = bxe[1:-1, 1, :]
        bxe[1:-1, -1, :] = bxe[1:-1, -2, :]

    if data.inp.mtd_vbnd[2] == 0:
        bxe[0, :, :] = bxe[-2, :, :]
        bxe[-1, :, :] = bxe[1, :, :]
    elif data.inp.mtd_vbnd[2] == 1:
        bxe[0, :, :] = -bxe[1, :, :]
        bxe[-1, :, :] = -bxe[-2, :, :]
    else:
        bxe[0, :, :] = bxe[1, :, :]
        bxe[-1, :, :] = bxe[-2, :, :]

    ebvalues[:, :, :, ielem] = 0.25 * (
        bxe[:-1, :-1, :] + bxe[1:, :-1, :] + bxe[:-1, 1:, :] + bxe[1:, 1:, :]
    )
    bxe = None

    # BY
    ielem = 4
    by = data.by[istep]
    # yz平面に1グリッド覆うように拡張する
    bye = np.empty(np.array(by.shape) + np.array([1, 0, 1]))
    bye[1:-1, :, 1:-1] = by[:-1, :, :-1]
    if data.inp.mtd_vbnd[2] == 0:
        bye[0, :, 1:-1] = bye[-2, :, 1:-1]
        bye[
            -1,
            :,
            1:-1,
        ] = bye[1, :, 1:-1]
    elif data.inp.mtd_vbnd[2] == 1:
        bye[0, :, 1:-1] = -bye[1, :, 1:-1]
        bye[-1, :, 1:-1] = -bye[-2, :, 1:-1]
    else:
        bye[0, :, 1:-1] = bye[1, :, 1:-1]
        bye[-1, :, 1:-1] = bye[-2, :, 1:-1]

    if data.inp.mtd_vbnd[0] == 0:
        bye[:, :, 0] = bye[:, :, -2]
        bye[:, :, -1] = bye[:, :, 1]
    elif data.inp.mtd_vbnd[0] == 1:
        bye[:, :, 0] = -bye[:, :, 1]
        bye[:, :, -1] = -bye[:, :, -2]
    else:
        bye[:, :, 0] = bye[:, :, 1]
        bye[:, :, -1] = bye[:, :, -2]

    ebvalues[:, :, :, ielem] = 0.25 * (
        bye[:-1, :, :-1] + bye[1:, :, :-1] + bye[:-1, :, 1:] + bye[1:, :, 1:]
    )
    bye = None

    # BZ
    ielem = 5
    bz = data.bz[istep]
    # xy平面に1グリッド覆うように拡張する
    bze = np.empty(np.array(bz.shape) + np.array([0, 1, 1]))
    bze[:, 1:-1, 1:-1] = bz[:, :-1, :-1]

    if data.inp.mtd_vbnd[0] == 0:
        bze[:, 1:-1, 0] = bze[:, 1:-1, -2]
        bze[:, 1:-1, -1] = bze[:, 1:-1, 1]
    elif data.inp.mtd_vbnd[0] == 1:
        bze[:, 1:-1, 0] = -bze[:, 1:-1, 1]
        bze[:, 1:-1, -1] = -bze[:, 1:-1, -2]
    else:
        bze[:, 1:-1, 0] = bze[:, 1:-1, 1]
        bze[:, 1:-1, -1] = bze[:, 1:-1, -2]

    if data.inp.mtd_vbnd[1] == 0:
        bze[:, 0, :] = bze[:, -2, :]
        bze[:, -1, :] = bze[:, 1, :]
    elif data.inp.mtd_vbnd[1] == 1:
        bze[:, 0, :] = -bze[:, 1, :]
        bze[:, -1, :] = -bze[:, -2, :]
    else:
        bze[:, 0, :] = bze[:, 1, :]
        bze[:, -1, :] = bze[:, -2, :]

    ebvalues[:, :, :, ielem] = 0.25 * (
        bze[:, :-1, :-1] + bze[:, 1:, :-1] + bze[:, :-1, 1:] + bze[:, 1:, 1:]
    )
    bze = None

    return ebvalues