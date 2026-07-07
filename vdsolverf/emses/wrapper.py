import os
import platform
from ctypes import *
from os import PathLike
from pathlib import Path
from typing import List, Literal, Tuple, Union

import emout
import numpy as np
from scipy.spatial.transform import Rotation

from ..core import Particle, VelocityRangeMap
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

_ELECTRIC_COMPONENTS = (
    ("ex", "rex", 2),
    ("ey", "rey", 1),
    ("ez", "rez", 0),
)


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
    tmp_input_suffix: Union[str, None] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

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
        tmp_input_suffix=tmp_input_suffix,
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
    tmp_input_suffix: Union[str, None] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
        tmp_input_suffix=tmp_input_suffix,
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
    n_threads: Union[int, None] = 1,
    tmp_input_suffix: Union[str, None] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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

    with TempolaryInput(data, suffix=tmp_input_suffix) as tmpinp:
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
        _n_threads = c_int(1 if n_threads is None else n_threads)

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
    tmp_input_suffix: Union[str, None] = None,
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
        tmp_input_suffix=tmp_input_suffix,
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
    tmp_input_suffix: Union[str, None] = None,
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
    with TempolaryInput(data, suffix=tmp_input_suffix) as tmpinp:
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


def estimate_velocity_range_map(
    directory: PathLike,
    ispec: int,
    istep: int,
    dt: float = 0.25,
    max_step: int = 1000,
    use_adaptive_dt: bool = False,
    coverage_sigma: Union[float, None] = 4.0,
    coverage_mode: Literal["sigma", "relative_density"] = "sigma",
    eps_rel: float = 1e-6,
    safety_factor: float = 1.25,
    max_probability_types: int = 100,
    source_samples_per_cell: int = 1,
    velocity_sample_mode: Literal["ellipsoid"] = "ellipsoid",
    minimum_count: int = 4,
    collect_moments: bool = False,
    show_progress: bool = True,
    accumulator_cache_size: int = 20000,
    system: Literal["auto", "linux", "darwin", "windows"] = "auto",
    library_path: PathLike = None,
    n_threads: Union[int, None] = None,
    tmp_input_suffix: Union[str, None] = None,
) -> VelocityRangeMap:
    """Estimate per-cell velocity ranges from EMSES source envelopes.

    By default this follows the legacy backtrace/probability APIs and uses the
    ordinary time-step interpretation. With ``use_adaptive_dt=True``, ``dt`` is
    treated as an approximate grid-distance limit per trace step, not as the
    EMSES simulation time step.
    """
    n_threads = n_threads or int(os.environ.get("OMP_NUM_THREADS", default="1"))

    dll = _load_dll(system, library_path)

    return estimate_velocity_range_map_dll(
        directory=directory,
        ispec=ispec,
        istep=istep,
        dt=dt,
        max_step=max_step,
        use_adaptive_dt=use_adaptive_dt,
        coverage_sigma=_resolve_coverage_sigma(
            coverage_sigma,
            coverage_mode=coverage_mode,
            eps_rel=eps_rel,
        ),
        safety_factor=safety_factor,
        max_probability_types=max_probability_types,
        source_samples_per_cell=source_samples_per_cell,
        velocity_sample_mode=velocity_sample_mode,
        minimum_count=minimum_count,
        collect_moments=collect_moments,
        show_progress=show_progress,
        accumulator_cache_size=accumulator_cache_size,
        dll=dll,
        n_threads=n_threads,
        tmp_input_suffix=tmp_input_suffix,
    )


def estimate_velocity_range_map_dll(
    directory: PathLike,
    ispec: int,
    istep: int,
    dt: float,
    max_step: int,
    use_adaptive_dt: bool,
    coverage_sigma: float,
    safety_factor: float,
    max_probability_types: int,
    dll: Union[CDLL, "WinDLL"],
    source_samples_per_cell: int = 1,
    velocity_sample_mode: Literal["ellipsoid"] = "ellipsoid",
    minimum_count: int = 4,
    collect_moments: bool = False,
    show_progress: bool = True,
    accumulator_cache_size: int = 20000,
    n_threads: int = 1,
    tmp_input_suffix: Union[str, None] = None,
) -> VelocityRangeMap:
    dll.estimate_velocity_range_map.argtypes = [
        c_char_p,  # inppath
        c_int,  # length
        c_int,  # lx
        c_int,  # ly
        c_int,  # lz
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),  # ebvalues
        c_int,  # ispec
        c_double,  # dt
        c_double,  # coverage_sigma
        c_double,  # safety_factor
        c_int,  # max_step
        c_int,  # use_adaptive_dt
        c_int,  # max_probability_types
        c_int,  # source_samples_per_cell
        c_int,  # velocity_sample_mode
        c_int,  # minimum_count
        c_int,  # collect_moments
        c_int,  # show_progress
        c_int,  # accumulator_cache_size
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # vx_min
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # vx_max
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # vy_min
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # vy_max
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # vz_min
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # vz_max
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=3),  # count
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # weight_sum
        c_void_p,  # mean_v
        c_void_p,  # cov_v
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=3),  # status
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=3),  # confidence
        POINTER(c_int),  # n_threads
    ]
    dll.estimate_velocity_range_map.restype = None

    data = emout.Emout(directory)

    ebvalues = create_relocated_ebvalues(data, istep, ispec=ispec)
    nz, ny, nx = data.inp.nz, data.inp.ny, data.inp.nx

    vx_min = np.empty((nz, ny, nx), dtype=np.float64)
    vx_max = np.empty((nz, ny, nx), dtype=np.float64)
    vy_min = np.empty((nz, ny, nx), dtype=np.float64)
    vy_max = np.empty((nz, ny, nx), dtype=np.float64)
    vz_min = np.empty((nz, ny, nx), dtype=np.float64)
    vz_max = np.empty((nz, ny, nx), dtype=np.float64)
    count = np.empty((nz, ny, nx), dtype=np.int32)
    weight_sum = np.empty((nz, ny, nx), dtype=np.float64)
    mean_v = np.empty((nz, ny, nx, 3), dtype=np.float64) if collect_moments else None
    cov_v = np.empty((nz, ny, nx, 3, 3), dtype=np.float64) if collect_moments else None
    status = np.empty((nz, ny, nx), dtype=np.int32)
    confidence = np.empty((nz, ny, nx), dtype=np.float64)

    with TempolaryInput(data, suffix=tmp_input_suffix) as tmpinp:
        inppath = tmpinp.tmppath
        inppath_str = str(inppath.resolve())

        _inppath = create_string_buffer(inppath_str.encode())
        _length = c_int(len(inppath_str))
        _nx = c_int(nx)
        _ny = c_int(ny)
        _nz = c_int(nz)
        _ispec = c_int(ispec + 1)
        _dt = c_double(dt)
        _coverage_sigma = c_double(coverage_sigma)
        _safety_factor = c_double(safety_factor)
        _max_step = c_int(max_step)
        _use_adaptive_dt = c_int(1 if use_adaptive_dt else 0)
        _max_probability_types = c_int(max_probability_types)
        _source_samples_per_cell = c_int(source_samples_per_cell)
        _velocity_sample_mode = c_int(_velocity_sample_mode_code(velocity_sample_mode))
        _minimum_count = c_int(minimum_count)
        _collect_moments = c_int(1 if collect_moments else 0)
        _show_progress = c_int(1 if show_progress else 0)
        _accumulator_cache_size = c_int(accumulator_cache_size)
        _mean_v = mean_v.ctypes.data_as(c_void_p) if mean_v is not None else c_void_p()
        _cov_v = cov_v.ctypes.data_as(c_void_p) if cov_v is not None else c_void_p()
        _n_threads = c_int(n_threads)

        dll.estimate_velocity_range_map(
            _inppath,
            _length,
            _nx,
            _ny,
            _nz,
            ebvalues,
            _ispec,
            _dt,
            _coverage_sigma,
            _safety_factor,
            _max_step,
            _use_adaptive_dt,
            _max_probability_types,
            _source_samples_per_cell,
            _velocity_sample_mode,
            _minimum_count,
            _collect_moments,
            _show_progress,
            _accumulator_cache_size,
            vx_min,
            vx_max,
            vy_min,
            vy_max,
            vz_min,
            vz_max,
            count,
            weight_sum,
            _mean_v,
            _cov_v,
            status,
            confidence,
            byref(_n_threads),
        )

    invalid = count <= 0
    for values in (vx_min, vx_max, vy_min, vy_max, vz_min, vz_max):
        values[invalid] = np.nan

    return VelocityRangeMap(
        x_edges=np.arange(nx + 1, dtype=np.float64),
        y_edges=np.arange(ny + 1, dtype=np.float64),
        z_edges=np.arange(nz + 1, dtype=np.float64),
        vx_min=vx_min,
        vx_max=vx_max,
        vy_min=vy_min,
        vy_max=vy_max,
        vz_min=vz_min,
        vz_max=vz_max,
        count=count,
        weight_sum=weight_sum,
        mean_v=mean_v,
        cov_v=cov_v,
        status=status,
        confidence=confidence,
        metadata={
            "directory": str(data.directory),
            "ispec": ispec,
            "istep": istep,
            "dt": dt,
            "max_step": max_step,
            "use_adaptive_dt": use_adaptive_dt,
            "coverage_sigma": coverage_sigma,
            "safety_factor": safety_factor,
            "source_samples_per_cell": source_samples_per_cell,
            "velocity_sample_mode": velocity_sample_mode,
            "minimum_count": minimum_count,
            "collect_moments": collect_moments,
            "show_progress": show_progress,
            "accumulator_cache_size": accumulator_cache_size,
            "n_threads": n_threads,
        },
        directory=data.directory,
    )


def validate_and_expand_velocity_range_map(
    range_map: VelocityRangeMap,
    directory: PathLike,
    ispec: int,
    istep: int,
    dt: float,
    max_step: int,
    coarse_bins: Tuple[int, int, int] = (8, 4, 8),
    edge_threshold: float = 1e-3,
    expand_factor: float = 1.5,
    max_iter: int = 2,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: Literal["auto", "linux", "darwin", "windows"] = "auto",
    library_path: PathLike = None,
    n_threads: Union[int, None] = None,
    tmp_input_suffix: Union[str, None] = None,
) -> VelocityRangeMap:
    for _ in range(max_iter):
        particles, index = range_map.create_particles(coarse_bins)
        if not particles:
            break

        probabilities, _ = get_probabilities(
            directory=directory,
            ispec=ispec,
            istep=istep,
            particles=particles,
            dt=dt,
            max_step=max_step,
            use_adaptive_dt=use_adaptive_dt,
            max_probability_types=max_probability_types,
            system=system,
            library_path=library_path,
            n_threads=n_threads,
            tmp_input_suffix=tmp_input_suffix,
        )
        prob_grid = index.reshape(probabilities)
        expand_mask = _edge_probability_mask(prob_grid, edge_threshold)
        if not expand_mask.any():
            break

        range_map.expand_cells(expand_mask, expand_factor)
        range_map.status[expand_mask] = 3

    return range_map


def _edge_probability_mask(prob_grid: np.ndarray, edge_threshold: float) -> np.ndarray:
    nz, ny, nx = prob_grid.shape[:3]
    ret = np.zeros((nz, ny, nx), dtype=bool)

    edge = np.zeros(prob_grid.shape[3:], dtype=bool)
    edge[0, :, :] = True
    edge[-1, :, :] = True
    edge[:, 0, :] = True
    edge[:, -1, :] = True
    edge[:, :, 0] = True
    edge[:, :, -1] = True

    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                cell = prob_grid[iz, iy, ix]
                if np.isnan(cell).all():
                    continue
                cell_max = np.nanmax(cell)
                if not np.isfinite(cell_max) or cell_max <= 0:
                    continue
                edge_max = np.nanmax(cell[edge])
                ret[iz, iy, ix] = np.isfinite(edge_max) and (
                    edge_max / cell_max > edge_threshold
                )

    return ret


def _resolve_coverage_sigma(
    coverage_sigma: Union[float, None],
    *,
    coverage_mode: Literal["sigma", "relative_density"],
    eps_rel: float,
) -> float:
    if coverage_mode == "sigma":
        if coverage_sigma is None:
            raise ValueError("coverage_sigma is required when coverage_mode='sigma'")
        return float(coverage_sigma)
    if coverage_mode == "relative_density":
        if eps_rel <= 0 or eps_rel >= 1:
            raise ValueError("eps_rel must be between 0 and 1")
        return float(np.sqrt(-2.0 * np.log(eps_rel)))
    raise ValueError(f"Unsupported coverage_mode: {coverage_mode}")


def _velocity_sample_mode_code(mode: Literal["ellipsoid"]) -> int:
    if mode == "ellipsoid":
        return 0
    raise ValueError(f"Unsupported velocity_sample_mode: {mode}")


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
        fill_separated_electric_fields(
            ebvalues, data, istep, phibk, field_substeps_per_particle_step(data)
        )
    else:
        fill_relocated_electric_field(ebvalues, data, istep)

    return ebvalues


def fill_separated_electric_fields(
    ebvalues: np.ndarray,
    data: emout.Emout,
    istep: int,
    accumulated_potential: np.ndarray,
    substeps_per_particle_step: float,
):
    expected_shape = (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1)

    for component, (field_name, relocated_name, axis) in enumerate(
        _ELECTRIC_COMPONENTS
    ):
        accumulated_e = ebvalues[:, :, :, component + 6]
        write_potential_difference(
            accumulated_potential,
            axis=axis,
            out=accumulated_e,
            scale=substeps_per_particle_step,
        )

        relocated_space_e = ebvalues[:, :, :, component]
        try:
            total_e = load_required_grid_step(data, field_name, istep)
        except AttributeError:
            relocated_total_e = load_required_grid_step(data, relocated_name, istep)
            validate_grid_shape(relocated_name, relocated_total_e, expected_shape)
            relocate_electric_component(
                accumulated_e,
                axis=axis,
                btype=electric_boundary_type(data, axis),
                out=relocated_space_e,
            )
            np.subtract(relocated_total_e, relocated_space_e, out=relocated_space_e)
            continue

        validate_grid_shape(field_name, total_e, expected_shape)
        relocate_electric_component_difference(
            total_e,
            accumulated_e,
            axis=axis,
            btype=electric_boundary_type(data, axis),
            out=relocated_space_e,
        )


def fill_relocated_electric_field(
    ebvalues: np.ndarray, data: emout.Emout, istep: int
):
    expected_shape = (data.inp.nz + 1, data.inp.ny + 1, data.inp.nx + 1)

    for component, (_, relocated_name, _) in enumerate(_ELECTRIC_COMPONENTS):
        relocated_e = load_required_grid_step(data, relocated_name, istep)
        validate_grid_shape(relocated_name, relocated_e, expected_shape)
        ebvalues[:, :, :, component] = relocated_e


def validate_grid_shape(name: str, grid: np.ndarray, expected_shape: Tuple[int, ...]):
    if grid.shape != expected_shape:
        raise ValueError(
            f"{name} shape mismatch: expected {expected_shape}, got {grid.shape}"
        )


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


def write_potential_difference(
    potential: np.ndarray, axis: int, out: np.ndarray, scale: float
):
    out[:, :, :] = 0.0
    if potential.shape[axis] <= 1:
        return

    left = [slice(None)] * 3
    left[axis] = slice(0, -1)
    right = [slice(None)] * 3
    right[axis] = slice(1, None)
    np.subtract(potential[tuple(left)], potential[tuple(right)], out=out[tuple(left)])

    last = [slice(None)] * 3
    last[axis] = -1
    previous = [slice(None)] * 3
    previous[axis] = -2
    out[tuple(last)] = out[tuple(previous)]

    if scale != 1.0:
        out[:, :, :] *= scale


def relocate_electric_component(
    component: np.ndarray,
    axis: int,
    btype: Literal["periodic", "dirichlet", "neumann"],
    out: Union[np.ndarray, None] = None,
) -> np.ndarray:
    relocated = np.zeros_like(component, dtype=np.float64) if out is None else out
    relocated[:, :, :] = 0.0

    if component.shape[axis] <= 1:
        return relocated

    middle = axis_slice(axis, slice(1, -1))
    lower = axis_slice(axis, slice(None, -2))
    upper = axis_slice(axis, slice(1, -1))
    np.add(component[lower], component[upper], out=relocated[middle])
    relocated[middle] *= 0.5

    first = axis_slice(axis, 0)
    last = axis_slice(axis, -1)

    if btype == "periodic":
        np.add(
            component[axis_slice(axis, -2)],
            component[axis_slice(axis, 1)],
            out=relocated[first],
        )
        relocated[first] *= 0.5
        relocated[last] = relocated[first]
    elif btype == "neumann":
        relocated[first] = 0.0
        relocated[last] = 0.0
    else:
        relocated[first] = component[axis_slice(axis, 1)]
        relocated[last] = component[axis_slice(axis, -2)]

    return relocated


def relocate_electric_component_difference(
    total: np.ndarray,
    accumulated: np.ndarray,
    axis: int,
    btype: Literal["periodic", "dirichlet", "neumann"],
    out: np.ndarray,
) -> np.ndarray:
    out[:, :, :] = 0.0

    if total.shape[axis] <= 1:
        return out

    middle = axis_slice(axis, slice(1, -1))
    lower = axis_slice(axis, slice(None, -2))
    upper = axis_slice(axis, slice(1, -1))
    np.subtract(total[lower], accumulated[lower], out=out[middle])
    out[middle] += total[upper]
    out[middle] -= accumulated[upper]
    out[middle] *= 0.5

    first = axis_slice(axis, 0)
    last = axis_slice(axis, -1)

    if btype == "periodic":
        np.subtract(
            total[axis_slice(axis, -2)],
            accumulated[axis_slice(axis, -2)],
            out=out[first],
        )
        out[first] += total[axis_slice(axis, 1)]
        out[first] -= accumulated[axis_slice(axis, 1)]
        out[first] *= 0.5
        out[last] = out[first]
    elif btype == "neumann":
        out[first] = 0.0
        out[last] = 0.0
    else:
        np.subtract(
            total[axis_slice(axis, 1)],
            accumulated[axis_slice(axis, 1)],
            out=out[first],
        )
        np.subtract(
            total[axis_slice(axis, -2)],
            accumulated[axis_slice(axis, -2)],
            out=out[last],
        )

    return out


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
