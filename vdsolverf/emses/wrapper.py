import os
import platform
from ctypes import *
from os import PathLike
from pathlib import Path
from typing import List, Literal, Tuple, Union

import emout
import numpy as np
from scipy.spatial.transform import Rotation

from ..core import Particle
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

_DEFAULT_LIBRARY_PATHS = {
    "linux": VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX,
    "darwin": VDIST_SOLVER_FORTRAN_LIBRARY_PATH_DARWIN,
    "windows": VDIST_SOLVER_FORTRAN_LIBRARY_PATH_WINDOWS,
}


def _load_dll(
    system: Literal["auto", "linux", "darwin", "windows"],
    library_path: Union[PathLike, None],
) -> Union[CDLL, "WinDLL"]:
    if system == "auto":
        system = platform.system().lower()

    if system not in _DEFAULT_LIBRARY_PATHS:
        raise RuntimeError(f"This platform is not supported: {system}")

    resolved_path = Path(library_path) if library_path is not None else _DEFAULT_LIBRARY_PATHS[system]

    if system == "windows":
        return WinDLL(str(resolved_path.resolve()))  # type: ignore[name-defined]
    return CDLL(str(resolved_path))


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

    dll = _load_dll(system, library_path)

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

    dll = _load_dll(system, library_path)

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
    dll.get_backtraces.restype = None

    data = emout.Emout(directory)

    ebvalues = create_relocated_ebvalues(data, istep)

    npcls = len(particles)

    max_output_steps = int((max_step - 1) / output_interval + 2)
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

    dll = _load_dll(system, library_path)

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
