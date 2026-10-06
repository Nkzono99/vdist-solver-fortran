"""Optional MPI wrappers for EMSES backtrace solvers.

The public functions in :mod:`vdsolverf.emses.wrapper` stay serial/OpenMP
compatible.  This module adds particle-parallel MPI entry points that split the
particle list across ranks and call the existing ctypes backend on each rank.
"""

from __future__ import annotations

import os
import pickle
import shlex
import subprocess
import sys
import tempfile
from os import PathLike
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple, Union

import numpy as np

from ..core import Particle
from . import wrapper as serial_wrapper


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
    system: str = "auto",
    library_path: Union[PathLike, None] = None,
    comm: Any = None,
    root: int = 0,
    return_on_all: bool = True,
    prepare_fields: bool = True,
    *,
    use_electric_field: bool = True,
    use_magnetic_field: bool = True,
):
    """Run a single-particle backtrace under MPI.

    Only ``root`` evaluates the particle because a single trajectory does not
    expose particle parallelism.  By default the result is broadcast so callers
    can use the same code on every rank.
    """
    comm = _get_comm(comm)
    rank = comm.Get_rank()

    if prepare_fields:
        _prepare_emout_inputs(
            comm, directory, istep, ispec, root,
            use_electric_field=use_electric_field,
            use_magnetic_field=use_magnetic_field,
        )

    result = None
    if rank == root:
        try:
            result = serial_wrapper.get_backtrace(
                directory=directory,
                ispec=ispec,
                istep=istep,
                particle=particle,
                dt=dt,
                max_step=max_step,
                output_interval=output_interval,
                use_adaptive_dt=use_adaptive_dt,
                max_probability_types=max_probability_types,
                system=system,
                library_path=library_path,
                tmp_input_suffix=_rank_tmp_input_suffix(comm),
                use_electric_field=use_electric_field,
                use_magnetic_field=use_magnetic_field,
            )
            error = None
        except Exception as exc:
            error = _error_info(rank, exc)
    else:
        error = None
    _raise_if_any_rank_failed(comm, error, root)

    return _return_result(comm, result, root, return_on_all)


def get_backtraces(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: Sequence[Particle],
    dt: float,
    max_step: int,
    output_interval: int = 1,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: str = "auto",
    library_path: Union[PathLike, None] = None,
    n_threads: Union[int, None] = None,
    comm: Any = None,
    root: int = 0,
    return_on_all: bool = True,
    prepare_fields: bool = True,
    *,
    use_electric_field: bool = True,
    use_magnetic_field: bool = True,
):
    """Run multi-particle backtraces by splitting particles across MPI ranks."""
    comm = _get_comm(comm)
    rank = comm.Get_rank()
    size = comm.Get_size()
    particles = list(particles)
    local_slice = _rank_slice(len(particles), size, rank)
    local_particles = particles[local_slice]

    if prepare_fields:
        _prepare_emout_inputs(
            comm, directory, istep, ispec, root,
            use_electric_field=use_electric_field,
            use_magnetic_field=use_magnetic_field,
        )

    try:
        if local_particles:
            local_result = serial_wrapper.get_backtraces(
                directory=directory,
                ispec=ispec,
                istep=istep,
                particles=local_particles,
                dt=dt,
                max_step=max_step,
                output_interval=output_interval,
                use_adaptive_dt=use_adaptive_dt,
                max_probability_types=max_probability_types,
                system=system,
                library_path=library_path,
                n_threads=n_threads,
                tmp_input_suffix=_rank_tmp_input_suffix(comm),
                use_electric_field=use_electric_field,
                use_magnetic_field=use_magnetic_field,
            )
        else:
            local_result = _empty_backtraces(max_step, output_interval)
        error = None
    except Exception as exc:
        local_result = _empty_backtraces(max_step, output_interval)
        error = _error_info(rank, exc)
    _raise_if_any_rank_failed(comm, error, root)

    gathered = comm.gather((local_slice.start, local_result), root=root)

    result = None
    if rank == root:
        result = _combine_backtrace_results(gathered, max_step, output_interval)

    return _return_result(comm, result, root, return_on_all)


def get_probabilities(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: Sequence[Particle],
    dt: float,
    max_step: int,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: str = "auto",
    library_path: Union[PathLike, None] = None,
    n_threads: Union[int, None] = None,
    comm: Any = None,
    root: int = 0,
    return_on_all: bool = True,
    prepare_fields: bool = True,
    *,
    use_electric_field: bool = True,
    use_magnetic_field: bool = True,
):
    """Run probability evaluation by splitting particles across MPI ranks."""
    comm = _get_comm(comm)
    rank = comm.Get_rank()
    size = comm.Get_size()
    particles = list(particles)
    local_slice = _rank_slice(len(particles), size, rank)
    local_particles = particles[local_slice]

    if prepare_fields:
        _prepare_emout_inputs(
            comm, directory, istep, ispec, root,
            use_electric_field=use_electric_field,
            use_magnetic_field=use_magnetic_field,
        )

    try:
        if local_particles:
            local_result = serial_wrapper.get_probabilities(
                directory=directory,
                ispec=ispec,
                istep=istep,
                particles=local_particles,
                dt=dt,
                max_step=max_step,
                use_adaptive_dt=use_adaptive_dt,
                max_probability_types=max_probability_types,
                system=system,
                library_path=library_path,
                n_threads=n_threads,
                tmp_input_suffix=_rank_tmp_input_suffix(comm),
                use_electric_field=use_electric_field,
                use_magnetic_field=use_magnetic_field,
            )
        else:
            local_result = (np.empty(0, dtype=np.float64), [])
        error = None
    except Exception as exc:
        local_result = (np.empty(0, dtype=np.float64), [])
        error = _error_info(rank, exc)
    _raise_if_any_rank_failed(comm, error, root)

    gathered = comm.gather((local_slice.start, local_result), root=root)

    result = None
    if rank == root:
        result = _combine_probability_results(gathered)

    return _return_result(comm, result, root, return_on_all)


def srun_get_backtrace(
    directory: PathLike,
    ispec: int,
    istep: int,
    particle: Particle,
    dt: float,
    max_step: int,
    output_interval: int = 1,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: str = "auto",
    library_path: Union[PathLike, None] = None,
    *,
    ntasks: Union[int, None] = None,
    launcher: Union[str, Sequence[str]] = "srun",
    launcher_args: Union[Sequence[str], None] = None,
    cpus_per_task: Union[int, None] = None,
    python_executable: Union[PathLike, None] = None,
    env: Union[Mapping[str, str], None] = None,
    timeout: Union[float, None] = None,
    tmpdir: Union[PathLike, None] = None,
    root: int = 0,
    prepare_fields: bool = True,
    use_electric_field: bool = True,
    use_magnetic_field: bool = True,
):
    """Launch an MPI worker with ``srun`` and return a single backtrace."""
    spec = _worker_spec(
        "get_backtrace",
        root,
        dict(
            directory=_resolved_path(directory),
            ispec=ispec,
            istep=istep,
            particle=particle,
            dt=dt,
            max_step=max_step,
            output_interval=output_interval,
            use_adaptive_dt=use_adaptive_dt,
            max_probability_types=max_probability_types,
            system=system,
            library_path=_optional_resolved_path(library_path),
            prepare_fields=prepare_fields,
            use_electric_field=use_electric_field,
            use_magnetic_field=use_magnetic_field,
        ),
    )
    return _run_srun_worker(
        spec,
        ntasks=ntasks,
        launcher=launcher,
        launcher_args=launcher_args,
        cpus_per_task=cpus_per_task,
        python_executable=python_executable,
        env=env,
        timeout=timeout,
        tmpdir=tmpdir,
    )


def srun_get_backtraces(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: Sequence[Particle],
    dt: float,
    max_step: int,
    output_interval: int = 1,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: str = "auto",
    library_path: Union[PathLike, None] = None,
    n_threads: Union[int, None] = None,
    *,
    ntasks: Union[int, None] = None,
    launcher: Union[str, Sequence[str]] = "srun",
    launcher_args: Union[Sequence[str], None] = None,
    cpus_per_task: Union[int, None] = None,
    python_executable: Union[PathLike, None] = None,
    env: Union[Mapping[str, str], None] = None,
    timeout: Union[float, None] = None,
    tmpdir: Union[PathLike, None] = None,
    root: int = 0,
    prepare_fields: bool = True,
    use_electric_field: bool = True,
    use_magnetic_field: bool = True,
):
    """Launch an MPI worker with ``srun`` and return multi-backtrace arrays."""
    spec = _worker_spec(
        "get_backtraces",
        root,
        dict(
            directory=_resolved_path(directory),
            ispec=ispec,
            istep=istep,
            particles=list(particles),
            dt=dt,
            max_step=max_step,
            output_interval=output_interval,
            use_adaptive_dt=use_adaptive_dt,
            max_probability_types=max_probability_types,
            system=system,
            library_path=_optional_resolved_path(library_path),
            n_threads=n_threads,
            prepare_fields=prepare_fields,
            use_electric_field=use_electric_field,
            use_magnetic_field=use_magnetic_field,
        ),
    )
    return _run_srun_worker(
        spec,
        ntasks=ntasks,
        launcher=launcher,
        launcher_args=launcher_args,
        cpus_per_task=cpus_per_task,
        python_executable=python_executable,
        env=env,
        timeout=timeout,
        tmpdir=tmpdir,
    )


def srun_get_probabilities(
    directory: PathLike,
    ispec: int,
    istep: int,
    particles: Sequence[Particle],
    dt: float,
    max_step: int,
    use_adaptive_dt: bool = False,
    max_probability_types: int = 100,
    system: str = "auto",
    library_path: Union[PathLike, None] = None,
    n_threads: Union[int, None] = None,
    *,
    ntasks: Union[int, None] = None,
    launcher: Union[str, Sequence[str]] = "srun",
    launcher_args: Union[Sequence[str], None] = None,
    cpus_per_task: Union[int, None] = None,
    python_executable: Union[PathLike, None] = None,
    env: Union[Mapping[str, str], None] = None,
    timeout: Union[float, None] = None,
    tmpdir: Union[PathLike, None] = None,
    root: int = 0,
    prepare_fields: bool = True,
    use_electric_field: bool = True,
    use_magnetic_field: bool = True,
):
    """Launch an MPI worker with ``srun`` and return probability results."""
    spec = _worker_spec(
        "get_probabilities",
        root,
        dict(
            directory=_resolved_path(directory),
            ispec=ispec,
            istep=istep,
            particles=list(particles),
            dt=dt,
            max_step=max_step,
            use_adaptive_dt=use_adaptive_dt,
            max_probability_types=max_probability_types,
            system=system,
            library_path=_optional_resolved_path(library_path),
            n_threads=n_threads,
            prepare_fields=prepare_fields,
            use_electric_field=use_electric_field,
            use_magnetic_field=use_magnetic_field,
        ),
    )
    return _run_srun_worker(
        spec,
        ntasks=ntasks,
        launcher=launcher,
        launcher_args=launcher_args,
        cpus_per_task=cpus_per_task,
        python_executable=python_executable,
        env=env,
        timeout=timeout,
        tmpdir=tmpdir,
    )


def worker_main(argv: Union[Sequence[str], None] = None) -> int:
    """Entry point for ``python -m vdsolverf.emses.mpi_worker``."""
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2:
        print("usage: python -m vdsolverf.emses.mpi_worker SPEC RESULT", file=sys.stderr)
        return 2

    spec_path = Path(args[0])
    result_path = Path(args[1])

    with spec_path.open("rb") as fp:
        spec = pickle.load(fp)

    comm = _get_comm(None)
    root = int(spec.get("root", 0))
    operation = spec["operation"]
    kwargs = dict(spec["kwargs"])
    kwargs["comm"] = comm
    kwargs["root"] = root
    kwargs["return_on_all"] = False

    if operation == "get_backtrace":
        result = get_backtrace(**kwargs)
    elif operation == "get_backtraces":
        result = get_backtraces(**kwargs)
    elif operation == "get_probabilities":
        result = get_probabilities(**kwargs)
    else:
        raise ValueError(f"Unknown MPI worker operation: {operation}")

    if comm.Get_rank() == root:
        result_path.parent.mkdir(parents=True, exist_ok=True)
        with result_path.open("wb") as fp:
            pickle.dump(result, fp, protocol=pickle.HIGHEST_PROTOCOL)

    return 0


def _load_mpi():
    try:
        from mpi4py import MPI
    except Exception as exc:
        raise RuntimeError(
            "The vdsolverf MPI backend requires mpi4py. "
            "Install it with `pip install vdist-solver-fortran[mpi]` "
            "or load an environment where mpi4py is available."
        ) from exc

    return MPI


def _get_comm(comm: Any):
    if comm is not None:
        return comm
    return _load_mpi().COMM_WORLD


def _rank_slice(length: int, size: int, rank: int) -> slice:
    start = length * rank // size
    end = length * (rank + 1) // size
    return slice(start, end)


def _rank_tmp_input_suffix(comm: Any) -> str:
    return f"mpi-r{comm.Get_rank()}-p{os.getpid()}"


def _prepare_emout_inputs(
    comm: Any,
    directory: PathLike,
    istep: int,
    ispec: Union[int, None],
    root: int,
    *,
    use_electric_field: bool = True,
    use_magnetic_field: bool = True,
):
    error = None
    if comm.Get_rank() == root:
        try:
            data = serial_wrapper.emout.Emout(directory)
            serial_wrapper.create_relocated_ebvalues(
                data, istep, ispec=ispec,
                use_electric_field=use_electric_field,
                use_magnetic_field=use_magnetic_field,
            )
        except Exception as exc:
            error = _error_info(root, exc)
    _raise_if_any_rank_failed(comm, error, root)


def _return_result(comm: Any, result: Any, root: int, return_on_all: bool):
    if return_on_all:
        return comm.bcast(result, root=root)
    if comm.Get_rank() == root:
        return result
    return None


def _error_info(rank: int, exc: Exception) -> Tuple[int, str, str]:
    return rank, type(exc).__name__, str(exc)


def _raise_if_any_rank_failed(
    comm: Any, local_error: Union[Tuple[int, str, str], None], root: int
):
    gathered = comm.gather(local_error, root=root)

    message = None
    if comm.Get_rank() == root:
        errors = [error for error in gathered if error is not None]
        if errors:
            message = "; ".join(
                f"rank {rank}: {exc_type}: {detail}"
                for rank, exc_type, detail in errors
            )

    message = comm.bcast(message, root=root)
    if message is not None:
        raise RuntimeError(f"vdsolverf MPI backend failed: {message}")


def _empty_backtraces(max_step: int, output_interval: int):
    max_output_steps = int((max_step - 1) / output_interval + 2)
    return (
        np.empty((0, max_output_steps), dtype=np.float64),
        np.empty(0, dtype=np.float64),
        np.empty((0, max_output_steps, 3), dtype=np.float64),
        np.empty((0, max_output_steps, 3), dtype=np.float64),
        np.empty(0, dtype=np.int32),
    )


def _combine_backtrace_results(
    gathered: Sequence[Tuple[int, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]],
    max_step: int,
    output_interval: int,
):
    ordered = [result for _, result in sorted(gathered, key=lambda item: item[0])]
    nonempty = [result for result in ordered if result[1].shape[0] > 0]
    if not nonempty:
        return _empty_backtraces(max_step, output_interval)

    return (
        np.concatenate([result[0] for result in nonempty], axis=0),
        np.concatenate([result[1] for result in nonempty], axis=0),
        np.concatenate([result[2] for result in nonempty], axis=0),
        np.concatenate([result[3] for result in nonempty], axis=0),
        np.concatenate([result[4] for result in nonempty], axis=0),
    )


def _combine_probability_results(
    gathered: Sequence[Tuple[int, Tuple[np.ndarray, List[Particle]]]]
):
    ordered = [result for _, result in sorted(gathered, key=lambda item: item[0])]
    probability_parts = [result[0] for result in ordered if result[0].size > 0]
    probabilities = (
        np.concatenate(probability_parts, axis=0)
        if probability_parts
        else np.empty(0, dtype=np.float64)
    )

    particles: List[Particle] = []
    for _, local_particles in ordered:
        particles.extend(local_particles)

    return probabilities, particles


def _worker_spec(operation: str, root: int, kwargs: Dict[str, Any]) -> Dict[str, Any]:
    return {"operation": operation, "root": root, "kwargs": kwargs}


def _run_srun_worker(
    spec: Dict[str, Any],
    *,
    ntasks: Union[int, None],
    launcher: Union[str, Sequence[str]],
    launcher_args: Union[Sequence[str], None],
    cpus_per_task: Union[int, None],
    python_executable: Union[PathLike, None],
    env: Union[Mapping[str, str], None],
    timeout: Union[float, None],
    tmpdir: Union[PathLike, None],
):
    ntasks = _resolve_ntasks(ntasks)
    executable = str(python_executable or sys.executable)

    with tempfile.TemporaryDirectory(dir=tmpdir) as tmp:
        tmp_path = Path(tmp)
        spec_path = tmp_path / "vdsolverf-mpi-spec.pkl"
        result_path = tmp_path / "vdsolverf-mpi-result.pkl"

        with spec_path.open("wb") as fp:
            pickle.dump(spec, fp, protocol=pickle.HIGHEST_PROTOCOL)

        command = _launcher_command(
            launcher=launcher,
            ntasks=ntasks,
            launcher_args=launcher_args,
            cpus_per_task=cpus_per_task,
            python_executable=executable,
            spec_path=spec_path,
            result_path=result_path,
        )
        run_env = _launcher_env(env, spec["kwargs"].get("n_threads"), cpus_per_task)

        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=run_env,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                "vdsolverf MPI launcher failed with exit code "
                f"{completed.returncode}\ncommand: {' '.join(command)}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        if not result_path.exists():
            raise RuntimeError(
                "vdsolverf MPI launcher finished without writing a result file: "
                f"{result_path}"
            )

        with result_path.open("rb") as fp:
            return pickle.load(fp)


def _resolve_ntasks(ntasks: Union[int, None]) -> int:
    if ntasks is None:
        ntasks = int(os.environ.get("VDSOLVERF_MPI_NTASKS", "1"))
    if ntasks < 1:
        raise ValueError("ntasks must be >= 1")
    return int(ntasks)


def _launcher_command(
    *,
    launcher: Union[str, Sequence[str]],
    ntasks: int,
    launcher_args: Union[Sequence[str], None],
    cpus_per_task: Union[int, None],
    python_executable: str,
    spec_path: Path,
    result_path: Path,
) -> List[str]:
    if isinstance(launcher, str):
        command = shlex.split(launcher)
    else:
        command = list(launcher)

    command.extend(["-n", str(ntasks)])
    if cpus_per_task is not None:
        command.append(f"--cpus-per-task={int(cpus_per_task)}")
    if launcher_args is not None:
        command.extend(str(arg) for arg in launcher_args)
    command.extend(
        [
            python_executable,
            "-m",
            "vdsolverf.emses.mpi_worker",
            str(spec_path),
            str(result_path),
        ]
    )
    return command


def _launcher_env(
    env: Union[Mapping[str, str], None],
    n_threads: Union[int, None],
    cpus_per_task: Union[int, None],
) -> Dict[str, str]:
    run_env = os.environ.copy()
    if cpus_per_task is not None:
        run_env.setdefault("OMP_NUM_THREADS", str(int(cpus_per_task)))
    if n_threads is not None:
        run_env["OMP_NUM_THREADS"] = str(int(n_threads))
    if env is not None:
        run_env.update({str(key): str(value) for key, value in env.items()})
    return run_env


def _resolved_path(path: PathLike) -> str:
    return str(Path(path).expanduser().resolve())


def _optional_resolved_path(path: Union[PathLike, None]) -> Union[str, None]:
    if path is None:
        return None
    return _resolved_path(path)
