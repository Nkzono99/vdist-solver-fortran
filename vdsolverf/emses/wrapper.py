import os
import platform
from ctypes import *
from os import PathLike
from pathlib import Path
from typing import List, Literal, Tuple, Union

import emout
import numpy as np
from scipy.spatial.transform import Rotation

from ..core import DustParticle, Particle
from .tmpolary_input import TempolaryInput

VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX = (
    Path(__file__).parent.parent / "libvdist-solver-fortran.so"
)

VDIST_SOLVER_FORTRAN_LIBRARY_PATH_DARWIN = (
    Path(__file__).parent.parent / "libvdist-solver-fortran.dylib"
)

VDIST_SOLVER_FORTRAN_LIBRARY_PATH_WINDOWS = (
    Path(__file__).parent.parent / "libvdist-solver-fortran.dll"
)


def get_backtrace(
    directory: PathLike,
    ispec: int,
    istep: int,
    particle: Particle,
    dt: float,
    max_step: int,
    output_interval: int = 1,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: Literal["auto", "linux", "darwin", "windows"] = "auto",
    library_path: PathLike = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    if system == "auto":
        system = platform.system().lower()

    if system == "linux":
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX
        dll = CDLL(library_path)
    elif system == "darwin":  # TODO: CDLLがこのプラットフォームで使えるのか要検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_DARWIN
        dll = CDLL(library_path)
    elif system == "windows":  # TODO: 実際に動作するのかは未検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_WINDOWS
        dll = WinDLL(library_path) # type:ignore
    else:
        raise RuntimeError(f"This platform is not supported: {system}")

    result = get_backtraces_dll(
        directory=directory,
        ispec=ispec,
        istep=istep,
        particles=[particle],
        dt=dt,
        max_step=max_step,
        output_interval=output_interval,
        use_adaptive_dt=use_adaptive_dt,
        max_probability_types=max_probability_types,
        dll=dll,
        n_threads=1,
    )

    ts, probabilities, positions_list, velocities_list, last_indexes = result

    # For some reason, it crashes when I try to close it.
    # handle = dll._handle

    # if os == "linux":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif os == "darwin":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif os == "windows":
    #     windll.kernel32.FreeLibrary(handle)

    return (
        ts[0, : last_indexes[0]],
        probabilities[0],
        positions_list[0, : last_indexes[0], :].copy(),
        velocities_list[0, : last_indexes[0], :].copy(),
    )


def get_backtraces(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: List[Particle],
    dt: float,
    max_step: int,
    output_interval: int = 1,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: Literal["auto", "linux", "darwin", "windows"] = "auto",
    library_path: PathLike = None,
    n_threads: Union[int, None] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_threads = n_threads or int(os.environ.get("OMP_NUM_THREADS", default="1"))

    if system == "auto":
        system = platform.system().lower()

    if system == "linux":
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX
        dll = CDLL(library_path)
    elif system == "darwin":  # TODO: CDLLがこのプラットフォームで使えるのか要検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_DARWIN
        dll = CDLL(library_path)
    elif system == "windows":  # TODO: 実際に動作するのかは未検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_WINDOWS
        dll = WinDLL(library_path) # type:ignore
    else:
        raise RuntimeError(f"This platform is not supported: {system}")

    result = get_backtraces_dll(
        directory=directory,
        ispec=ispec,
        istep=istep,
        particles=particles,
        dt=dt,
        max_step=max_step,
        output_interval=output_interval,
        use_adaptive_dt=use_adaptive_dt,
        max_probability_types=max_probability_types,
        dll=dll,
        n_threads=n_threads,
    )

    # For some reason, it crashes when I try to close it.
    # handle = dll._handle

    # if os == "linux":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif os == "darwin":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif os == "windows":
    #     windll.kernel32.FreeLibrary(handle)

    return result


def get_backtraces_dll(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: List[Particle],
    dt: float,
    max_step: int,
    output_interval: int,
    use_adaptive_dt: bool,
    max_probability_types: int,
    dll: Union[CDLL, "WinDLL"],
    n_threads: Union[int, None] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    dll.get_backtraces.argtypes = [
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
        c_double,  # dt
        c_int,  # max_step
        c_int,  # output_interval
        c_int,  # use_adaptive_dt
        c_int,  # max_probability_types
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # return_ts
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),  # return_probability
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # return_positions
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # return_velocities
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1),  # return_last_step
        POINTER(c_int),  # n_threads
    ]
    dll.get_probabilities.restype = None

    data = emout.Emout(directory)

    ebvalues = create_relocated_ebvalues(data, istep)

    npcls = len(particles)

    max_output_steps = int((max_step - 1)/output_interval + 2)
    return_ts = np.empty((npcls, max_output_steps), dtype=np.float64)
    return_probabilities = np.empty(npcls, dtype=np.float64)
    return_positions = np.empty((npcls, max_output_steps, 3), dtype=np.float64)
    return_velocities = np.empty((npcls, max_output_steps, 3), dtype=np.float64)
    return_last_indexes = np.empty(npcls, dtype=np.int32)

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
        _npcls = c_int(npcls)
        _dt = c_double(dt)
        _max_step = c_int(max_step)
        _output_interval = c_int(output_interval)
        _use_adaptive_dt = c_int(1 if use_adaptive_dt else 0)
        _max_probability_types = c_int(max_probability_types)
        _n_threads = c_int(n_threads)

        dll.get_backtraces(
            _inppath,
            _length,
            _nx,
            _ny,
            _nz,
            ebvalues,
            _ispec,
            _npcls,
            positions,
            velocities,
            _dt,
            _max_step,
            _output_interval,
            _use_adaptive_dt,
            _max_probability_types,
            return_ts,
            return_probabilities,
            return_positions,
            return_velocities,
            return_last_indexes,
            byref(_n_threads),
        )

    return_probabilities[return_probabilities == -1] = np.nan

    return (
        return_ts,
        return_probabilities,
        return_positions,
        return_velocities,
        return_last_indexes,
    )


def get_probabilities(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: List[Particle],
    dt: float,
    max_step: int,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: Literal["auto", "linux", "darwin", "windows"] = "auto",
    library_path: PathLike = None,
    n_threads: Union[int, None] = None,
) -> Tuple[np.ndarray, List[Particle]]:
    n_threads = n_threads or int(os.environ.get("OMP_NUM_THREADS", default="1"))

    if system == "auto":
        system = platform.system().lower()

    if system == "linux":
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX
        dll = CDLL(library_path)
    elif system == "darwin":  # TODO: CDLLがこのプラットフォームで使えるのか要検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_DARWIN
        dll = CDLL(library_path)
    elif system == "windows":  # TODO: 実際に動作するのかは未検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_WINDOWS
        dll = WinDLL(library_path) # type:ignore
    else:
        raise RuntimeError(f"This platform is not supported: {system}")

    result = get_probabilities_dll(
        directory=directory,
        ispec=ispec,
        istep=istep,
        particles=particles,
        dt=dt,
        max_step=max_step,
        use_adaptive_dt=use_adaptive_dt,
        max_probability_types=max_probability_types,
        dll=dll,
        n_threads=n_threads,
    )

    # For some reason, it crashes when I try to close it.
    # handle = dll._handle

    # if system == "linux":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif system == "darwin":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif system == "windows":
    #     windll.kernel32.FreeLibrary(handle)

    return result


def get_probabilities_dll(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: List[Particle],
    dt: float,
    max_step: int,
    use_adaptive_dt: bool,
    max_probability_types: int,
    dll: Union[CDLL, "WinDLL"],
    n_threads: int = 1,
) -> Tuple[np.ndarray, List[Particle]]:
    dll.get_probabilities.argtypes = [
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
        c_double,  # dt
        c_int,  # max_step
        c_int,  # use_adaptive_dt
        c_int,  # max_probability_types
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),  # return_probabilities
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # return_positions
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # return_velocities
        POINTER(c_int),  # n_threads
    ]
    dll.get_probabilities.restype = None

    data = emout.Emout(directory)

    ebvalues = create_relocated_ebvalues(data, istep)

    npcls = len(particles)
    return_probabilities = np.empty(npcls, dtype=np.float64)
    return_positions = np.empty((npcls, 3), dtype=np.float64)
    return_velocities = np.empty((npcls, 3), dtype=np.float64)

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
        _nparticles = c_int(npcls)
        _dt = c_double(dt)
        _max_step = c_int(max_step)
        _use_adaptive_dt = c_int(1 if use_adaptive_dt else 0)
        _max_probability_types = c_int(max_probability_types)
        _n_threads = c_int(n_threads)

        dll.get_probabilities(
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
            _dt,
            _max_step,
            _use_adaptive_dt,
            _max_probability_types,
            return_probabilities,
            return_positions,
            return_velocities,
            byref(_n_threads),
        )

    return_particles = [
        Particle(pos, vel) for pos, vel in zip(return_positions, return_velocities)
    ]

    return_probabilities[return_probabilities == -1] = np.nan

    return return_probabilities, return_particles


def get_dust_backtrace(
    directory: PathLike,
    istep: int,
    dust: DustParticle,
    dt: float,
    max_step: int,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    os: Literal["auto", "linux", "darwin", "windows"] = "auto",
    library_path: PathLike = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if os == "auto":
        os = platform.system().lower()

    if os == "linux":
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX
        dll = CDLL(library_path)
    elif os == "darwin":  # TODO: CDLLがこのプラットフォームで使えるのか要検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_DARWIN
        dll = CDLL(library_path)
    elif os == "windows":  # TODO: 実際に動作するのかは未検証
        library_path = library_path or VDIST_SOLVER_FORTRAN_LIBRARY_PATH_WINDOWS
        dll = WinDLL(library_path) # type:ignore
    else:
        raise RuntimeError(f"This platform is not supported: {os}")

    result = get_dust_backtrace_dll(
        directory=directory,
        istep=istep,
        dust=dust,
        dt=dt,
        max_step=max_step,
        use_adaptive_dt=use_adaptive_dt,
        max_probability_types=max_probability_types,
        dll=dll,
    )

    # For some reason, it crashes when I try to close it.
    # handle = dll._handle

    # if os == "linux":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif os == "darwin":
    #     cdll.LoadLibrary("libdl.so").dlclose(handle)
    # elif os == "windows":
    #     windll.kernel32.FreeLibrary(handle)

    return result


def get_dust_backtrace_dll(
    directory: PathLike,
    istep: int,
    dust: DustParticle,
    dt: float,
    max_step: int,
    use_adaptive_dt: bool,
    max_probability_types: int,
    dll: Union[CDLL, "WinDLL"],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dll.get_backtrace_dust.argtypes = [
        c_char_p,  # inppath
        c_int,  # length
        c_int,  # lx
        c_int,  # ly
        c_int,  # lz
        c_int,  # nspec
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),  # ebvalues
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),  # current_values
        c_double,  # charge
        c_double,  # mass
        c_double,  # radius
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),  # positions
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),  # velocities
        c_double,  # dt
        c_int,  # max_step
        c_int,  # use_adaptive_dt
        c_int,  # max_probability_types
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),  # return_ts
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),  # return_charges
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # return_positions
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2),  # return_velocities
        POINTER(c_int),  # return_last_step
    ]
    dll.get_probabilities.restype = None

    data = emout.Emout(directory)

    ebvalues = create_relocated_ebvalues(data, istep)
    current_values = create_relocated_current_values(data, istep)

    return_ts = np.empty(max_step, dtype=np.float64)
    return_charges = np.empty(max_step, dtype=np.float64)
    return_positions = np.empty((max_step, 3), dtype=np.float64)
    return_velocities = np.empty((max_step, 3), dtype=np.float64)

    position = np.array(dust.pos, dtype=np.float64)
    velocity = np.array(dust.vel, dtype=np.float64)

    with TempolaryInput(data) as tmpinp:
        inppath = tmpinp.tmppath
        inppath_str = str(inppath.resolve())

        _inppath = create_string_buffer(inppath_str.encode())
        _length = c_int(len(inppath_str))
        _nx = c_int(data.inp.nx)
        _ny = c_int(data.inp.ny)
        _nz = c_int(data.inp.nz)
        _nspec = c_int(data.inp.nspec)
        _charge = c_double(dust.charge)
        _mass = c_double(dust.mass)
        _radius = c_double(dust.radius)
        _dt = c_double(dt)
        _max_step = c_int(max_step)
        _use_adaptive_dt = c_int(1 if use_adaptive_dt else 0)
        _max_probability_types = c_int(max_probability_types)
        _return_last_index = c_int()

        dll.get_backtrace_dust(
            _inppath,
            _length,
            _nx,
            _ny,
            _nz,
            _nspec,
            ebvalues,
            current_values,
            _charge,
            _mass,
            _radius,
            position,
            velocity,
            _dt,
            _max_step,
            _use_adaptive_dt,
            _max_probability_types,
            return_ts,
            return_charges,
            return_positions,
            return_velocities,
            byref(_return_last_index),
        )
        return_last_index = _return_last_index.value

    ts = return_ts[:return_last_index].copy()
    charges = return_charges[:return_last_index].copy()
    positions = return_positions[:return_last_index].copy()
    velocities = return_velocities[:return_last_index].copy()

    return ts, charges, positions, velocities


def create_relocated_ebvalues(data: emout.Emout, istep: int) -> np.ndarray:
    ebvalues = np.zeros(
        (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1, 6), dtype=np.float64
    )

    ebvalues[:, :, :, 0] = data.rex[istep, :, :, :]
    ebvalues[:, :, :, 1] = data.rey[istep, :, :, :]
    ebvalues[:, :, :, 2] = data.rez[istep, :, :, :]
    ebvalues[:, :, :, 3] = data.rbx[istep, :, :, :]
    ebvalues[:, :, :, 4] = data.rby[istep, :, :, :]
    ebvalues[:, :, :, 5] = data.rbz[istep, :, :, :]

    b0x, b0y, b0z = background_magnetic_field(data)

    ebvalues[:, :, :, 3] += b0x
    ebvalues[:, :, :, 4] += b0y
    ebvalues[:, :, :, 5] += b0z

    return ebvalues


def background_magnetic_field(data: emout.Emout) -> np.ndarray:
    if "wc" not in data.inp:
        return np.zeros(3)

    b0 = data.inp.wc / data.inp.qm[0]

    return rotate(np.array([0.0, 0.0, b0]), data.inp.phiz, data.inp.phixy)


def rotate(vec: np.ndarray, phiz_deg: float, phixy_deg: float) -> np.ndarray:
    rot = Rotation.from_euler("yz", [phiz_deg, phixy_deg], degrees=True)

    return rot.apply(vec)


def create_relocated_current_values(data: emout.Emout, istep: int) -> np.ndarray:
    current_values = np.zeros(
        (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1, 3 * data.inp.nspec),
        dtype=np.float64,
    )

    for ispec in range(data.inp.nspec):
        # EX
        ielem = 0
        jx = getattr(data, f"j{ispec+1}x")[istep]
        current_values[:, :, 1:-1, ispec * 3 + ielem] = 0.5 * (
            jx[:, :, :-2] + jx[:, :, 1:-1]
        )
        if data.inp.mtd_vbnd[0] in [0, 2]:
            current_values[:, :, 0, ispec * 3 + ielem] = 0
            current_values[:, :, -1, ispec * 3 + ielem] = 0
        else:
            current_values[:, :, 0, ispec * 3 + ielem] = jx[:, :, 0]
            current_values[:, :, -1, ispec * 3 + ielem] = jx[:, :, -1]

        # EY
        ielem = 1
        jy = getattr(data, f"j{ispec+1}y")[istep]
        current_values[:, 1:-1, :, ispec * 3 + ielem] = 0.5 * (
            jy[:, :-2, :] + jy[:, 1:-1, :]
        )
        if data.inp.mtd_vbnd[0] in [0, 2]:
            current_values[:, 0, :, ispec * 3 + ielem] = 0
            current_values[:, -1, :, ispec * 3 + ielem] = 0
        else:
            current_values[:, 0, :, ispec * 3 + ielem] = jy[:, 0, :]
            current_values[:, -1, :, ispec * 3 + ielem] = jy[:, -1, :]

        # EZ
        ielem = 2
        jz = getattr(data, f"j{ispec+1}z")[istep]
        current_values[1:-1, :, :, ispec * 3 + ielem] = 0.5 * (
            jz[:-2, :, :] + jz[1:-1, :, :]
        )
        if data.inp.mtd_vbnd[0] in [0, 2]:
            current_values[0, :, :, ispec * 3 + ielem] = 0
            current_values[-1, :, :, ispec * 3 + ielem] = 0
        else:
            current_values[0, :, :, ispec * 3 + ielem] = jz[0, :, :]
            current_values[-1, :, :, ispec * 3 + ielem] = jz[-1, :, :]

    return current_values
