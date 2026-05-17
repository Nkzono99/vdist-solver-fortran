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
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),  # ebvalues (9 components)
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

    ebvalues = create_relocated_ebvalues(data, istep, ispec=ispec)

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
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),  # ebvalues (9 components)
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

    ebvalues = create_relocated_ebvalues(data, istep, ispec=ispec)

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


def create_relocated_ebvalues(
    data: emout.Emout, istep: int, ispec: Union[int, None] = None
) -> np.ndarray:
    ebvalues = np.zeros(
        (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1, 9), dtype=np.float64
    )

    ebvalues[:, :, :, 3] = data.rbx[istep, :, :, :]
    ebvalues[:, :, :, 4] = data.rby[istep, :, :, :]
    ebvalues[:, :, :, 5] = data.rbz[istep, :, :, :]

    b0x, b0y, b0z = background_magnetic_field(data)

    ebvalues[:, :, :, 3] += b0x
    ebvalues[:, :, :, 4] += b0y
    ebvalues[:, :, :, 5] += b0z

    phibk = load_accumulated_potential(data, istep, ispec)
    if phibk is not None:
        accumulated_e = accumulated_electric_field_from_potential(
            phibk, field_substeps_per_particle_step(data)
        )
        ebvalues[:, :, :, 0:3] = create_relocated_space_electric_field(
            data, istep, accumulated_e
        )
        ebvalues[:, :, :, 6:9] = accumulated_e
    else:
        ebvalues[:, :, :, 0:3] = load_relocated_electric_field(data, istep)

    return ebvalues


def create_relocated_space_electric_field(
    data: emout.Emout, istep: int, accumulated_electric_field: np.ndarray
) -> np.ndarray:
    try:
        space_electric_field = (
            load_electric_field(data, istep) - accumulated_electric_field
        )
    except AttributeError:
        relocated_total_e = load_relocated_electric_field(data, istep)
        relocated_accumulated_e = relocate_electric_field(data, accumulated_electric_field)
        return relocated_total_e - relocated_accumulated_e

    return relocate_electric_field(data, space_electric_field)


def load_electric_field(data: emout.Emout, istep: int) -> np.ndarray:
    expected_shape = (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1)
    electric_field = np.zeros(expected_shape + (3,), dtype=np.float64)

    for component, name in enumerate(("ex", "ey", "ez")):
        grid = load_required_grid_step(data, name, istep)
        if grid.shape != expected_shape:
            raise ValueError(
                f"{name} shape mismatch: expected {expected_shape}, got {grid.shape}"
            )
        electric_field[:, :, :, component] = grid

    return electric_field


def load_relocated_electric_field(data: emout.Emout, istep: int) -> np.ndarray:
    expected_shape = (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1)
    electric_field = np.zeros(expected_shape + (3,), dtype=np.float64)

    for component, name in enumerate(("rex", "rey", "rez")):
        grid = load_required_grid_step(data, name, istep)
        if grid.shape != expected_shape:
            raise ValueError(
                f"{name} shape mismatch: expected {expected_shape}, got {grid.shape}"
            )
        electric_field[:, :, :, component] = grid

    return electric_field


def load_accumulated_potential(
    data: emout.Emout, istep: int, ispec: Union[int, None]
) -> Union[np.ndarray, None]:
    for name in accumulated_potential_names(ispec):
        potential = load_optional_grid_step(data, name, istep)
        if potential is None:
            continue

        expected_shape = (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1)
        if potential.shape != expected_shape:
            raise ValueError(
                f"{name} shape mismatch: expected {expected_shape}, got {potential.shape}"
            )
        return potential

    return None


def accumulated_potential_names(ispec: Union[int, None]) -> List[str]:
    names = []

    if ispec is not None:
        species = ispec + 1
        names.extend([f"phibksp{species}", f"phibksp{species:02d}"])

    names.extend(["phibksp", "phibk"])

    return list(dict.fromkeys(names))


def load_optional_grid_step(
    data: emout.Emout, name: str, istep: int
) -> Union[np.ndarray, None]:
    try:
        series = getattr(data, name)
    except AttributeError:
        return None

    array = np.asarray(series[istep, :, :, :], dtype=np.float64)
    if array.ndim != 3:
        raise ValueError(f"{name} must be a 3D grid at one step, got ndim={array.ndim}")

    return array


def load_required_grid_step(data: emout.Emout, name: str, istep: int) -> np.ndarray:
    array = load_optional_grid_step(data, name, istep)
    if array is None:
        raise AttributeError(f"Required EMSES field '{name}' was not found")

    return array


def accumulated_electric_field_from_potential(
    potential: np.ndarray, substeps_per_particle_step: float
) -> np.ndarray:
    electric_field = np.zeros(potential.shape + (3,), dtype=np.float64)

    electric_field[:, :, :, 0] = potential_difference(potential, axis=2)
    electric_field[:, :, :, 1] = potential_difference(potential, axis=1)
    electric_field[:, :, :, 2] = potential_difference(potential, axis=0)

    electric_field *= substeps_per_particle_step

    return electric_field


def potential_difference(potential: np.ndarray, axis: int) -> np.ndarray:
    difference = np.zeros_like(potential, dtype=np.float64)

    if potential.shape[axis] <= 1:
        return difference

    left = [slice(None)] * 3
    left[axis] = slice(0, -1)
    right = [slice(None)] * 3
    right[axis] = slice(1, None)
    difference[tuple(left)] = potential[tuple(left)] - potential[tuple(right)]

    last = [slice(None)] * 3
    last[axis] = -1
    previous = [slice(None)] * 3
    previous[axis] = -2
    difference[tuple(last)] = difference[tuple(previous)]

    return difference


def relocate_electric_field(
    data: emout.Emout, electric_field: np.ndarray
) -> np.ndarray:
    relocated = np.zeros_like(electric_field, dtype=np.float64)

    for component, axis in enumerate((2, 1, 0)):
        relocated[:, :, :, component] = relocate_electric_component(
            electric_field[:, :, :, component],
            axis=axis,
            btype=electric_boundary_type(data, axis),
        )

    return relocated


def relocate_electric_component(
    component: np.ndarray, axis: int, btype: Literal["periodic", "dirichlet", "neumann"]
) -> np.ndarray:
    relocated = np.zeros_like(component, dtype=np.float64)

    if component.shape[axis] <= 1:
        return relocated

    middle = axis_slice(axis, slice(1, -1))
    lower = axis_slice(axis, slice(None, -2))
    upper = axis_slice(axis, slice(1, -1))
    relocated[middle] = 0.5 * (component[lower] + component[upper])

    first = axis_slice(axis, 0)
    last = axis_slice(axis, -1)

    if btype == "periodic":
        periodic_value = 0.5 * (
            component[axis_slice(axis, -2)] + component[axis_slice(axis, 1)]
        )
        relocated[first] = periodic_value
        relocated[last] = periodic_value
    elif btype == "neumann":
        relocated[first] = 0.0
        relocated[last] = 0.0
    else:
        relocated[first] = component[axis_slice(axis, 1)]
        relocated[last] = component[axis_slice(axis, -2)]

    return relocated


def axis_slice(axis: int, axis_index) -> Tuple:
    slices = [slice(None)] * 3
    slices[axis] = axis_index

    return tuple(slices)


def electric_boundary_type(
    data: emout.Emout, axis: int
) -> Literal["periodic", "dirichlet", "neumann"]:
    btypes = ("periodic", "dirichlet", "neumann")

    try:
        boundary_code = int(data.inp.mtd_vbnd[2 - axis])
        return btypes[boundary_code]
    except (AttributeError, IndexError, TypeError, ValueError):
        return "dirichlet"


def field_substeps_per_particle_step(data: emout.Emout) -> float:
    for name in ("field_substeps_per_particle_step", "mltstp"):
        if name in data.inp:
            return float(getattr(data.inp, name))

    return 1.0


def background_magnetic_field(data: emout.Emout) -> np.ndarray:
    if "wc" not in data.inp:
        return np.zeros(3)

    b0 = data.inp.wc / data.inp.qm[0]

    return rotate(np.array([0.0, 0.0, b0]), data.inp.phiz, data.inp.phixy)


def rotate(vec: np.ndarray, phiz_deg: float, phixy_deg: float) -> np.ndarray:
    rot = Rotation.from_euler("yz", [phiz_deg, phixy_deg], degrees=True)

    return rot.apply(vec)
