"""Prepare standalone current batches and collect completed chunk arrays."""
import json
import os
from pathlib import Path

import f90nml
import numpy as np

from ..core.dust_current import DustCurrentResult
from ..core.velocity_range_map import VelocityRangeMap
from .field_bundle import FieldBundle, read_field_bundle


_DEFAULTS = dict(
    mode="dust_currents",
    max_step=30_000, use_adaptive_dt=False, max_probability_types=100,
    scout_bins=[9, 9, 9], max_depth=6, max_samples_per_cell=50_000,
    max_leaves_per_cell=4096, refine_threshold_rel=1e-4, edge_threshold_rel=1e-4,
    expand_factor=1.5, max_expansions=3, chunk_size=16, output_dir="out",
    use_electric_field=True, use_magnetic_field=True,
)
_INTEGER_MIN = dict(ispec=0, max_step=1, max_probability_types=1, max_depth=0,
                    max_samples_per_cell=1, max_leaves_per_cell=1, max_expansions=0,
                    chunk_size=1, n_threads=1)
_BOOLEANS = ("use_adaptive_dt", "use_electric_field", "use_magnetic_field")


def _validate(config):
    unknown = set(config) - set(_DEFAULTS) - {"dt", "ispec", "bundle_dir", "n_threads"}
    if unknown:
        raise ValueError(f"unknown batch settings: {sorted(unknown)}")
    if config["mode"] != "dust_currents":
        raise ValueError("expected mode='dust_currents'")
    for name, minimum in _INTEGER_MIN.items():
        if name not in config:
            continue
        value = config[name]
        if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
                or not minimum <= value < np.iinfo(np.int32).max):
            raise ValueError(f"{name} must be an integer >= {minimum} within the int32 range")
        config[name] = int(value)
    bins = np.asarray(config["scout_bins"])
    if bins.shape != (3,) or bins.dtype.kind not in "iu" or np.any(bins < 2):
        raise ValueError("scout_bins must contain three integers >= 2")
    if np.prod(bins.astype(object)) > np.iinfo(np.int32).max:
        raise ValueError("scout_bins product exceeds the int32 range")
    config["scout_bins"] = [int(b) for b in bins]
    for name, minimum in (("dt", 0), ("refine_threshold_rel", 0),
                          ("edge_threshold_rel", 0), ("expand_factor", 1)):
        value = float(config[name])
        if not np.isfinite(value) or value < minimum or (name == "dt" and value == 0):
            raise ValueError(f"invalid batch {name}")
        config[name] = value
    for name in _BOOLEANS:
        if not isinstance(config[name], (bool, np.bool_)):
            raise ValueError(f"{name} must be a boolean")
        config[name] = bool(config[name])
    for name in ("bundle_dir", "output_dir"):
        value = os.fspath(config[name])
        if not value or len(value.encode()) > 4095 or any(c in value for c in "\x00\n\r"):
            raise ValueError(f"invalid batch {name}")
        config[name] = value


def _batch_path(path):
    path = Path(path).resolve()
    return path / "batch.nml" if path.is_dir() else path


def _read_inputs(path):
    nml_path = _batch_path(path)
    config = dict(_DEFAULTS)
    config.update(f90nml.read(nml_path)["vdsolverf_batch"])
    if "dt" not in config or "ispec" not in config or "bundle_dir" not in config:
        raise ValueError("batch.nml requires bundle_dir, ispec and dt")
    _validate(config)
    targets = _load_npy(nml_path.parent / "targets.npy", "<f8")
    psi = _load_npy(nml_path.parent / "psi.npy", "<f8")
    _validate_arrays(targets, psi)
    return nml_path, config, targets, psi


def _validate_arrays(targets, psi):
    if targets.ndim != 2 or targets.shape[1] != 9 or not len(targets):
        raise ValueError("targets must have nonempty shape (n, 9)")
    if not np.all(np.isfinite(targets[:, :3])):
        raise ValueError("target positions must be finite")
    # Let the solver report status 5 for invalid/nonfinite velocity boxes.
    if psi.ndim != 1 or not len(psi) or not np.all(np.isfinite(psi)):
        raise ValueError("psi must be a nonempty finite one-dimensional array")
    if len(targets) > np.iinfo(np.int32).max or len(psi) > np.iinfo(np.int32).max:
        raise ValueError("targets/psi exceed the native int32 index range")


def _save_npy(path, values):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.lib.format.write_array(stream, values, version=(1, 0), allow_pickle=False)
    temporary.replace(path)


def _load_npy(path, dtype, expected=None):
    with path.open("rb") as stream:
        version = np.lib.format.read_magic(stream)
        if version == (1, 0):
            shape, fortran_order, actual_dtype = np.lib.format.read_array_header_1_0(stream)
        elif version == (2, 0):
            shape, fortran_order, actual_dtype = np.lib.format.read_array_header_2_0(stream)
        else:
            raise ValueError(f"unsupported npy version in {path}")
        if fortran_order or actual_dtype.str != dtype or (expected is not None and shape != expected):
            raise ValueError(f"invalid chunk/input {path}: expected C-order {dtype} {expected}")
    values = np.load(path, allow_pickle=False)
    return values


def write_batch(batch_dir, bundle_dir, ispec, targets, psi, **solver):
    """Write batch.nml, targets.npy and psi.npy; return the batch.nml Path.

    targets is (n,9): x,y,z,vxmin,vxmax,vymin,vymax,vzmin,vzmax, in EMSES
    units. ispec is zero-based. dt defaults to the bundle's namelist dt;
    other settings follow get_dust_currents, with chunk_size=16. Omit
    n_threads to use OMP_NUM_THREADS at execution time. Paths passed here
    resolve from the caller's working directory and are stored relative to
    batch.nml. Keep inputs immutable when resuming completed chunks.
    """
    batch_dir = Path(batch_dir).resolve()
    bundle = bundle_dir if isinstance(bundle_dir, FieldBundle) else read_field_bundle(bundle_dir)
    bundle_path = bundle.namelist_path.parent.resolve()
    config = dict(_DEFAULTS)
    config.update(solver)
    if config.get("n_threads", 1) is None:
        config.pop("n_threads")
    config["ispec"] = ispec
    config["bundle_dir"] = os.path.relpath(bundle_path, batch_dir)
    if "dt" not in config:
        config["dt"] = f90nml.read(bundle.namelist_path)["tmgrid"]["dt"]
    _validate(config)
    if config["ispec"] not in bundle.metadata["ispecs"]:
        raise ValueError(f"species {ispec} is not present in the field bundle")
    for name in ("use_electric_field", "use_magnetic_field"):
        if config[name] != bundle.metadata[name]:
            raise ValueError(f"{name} must match the field bundle")
    targets = np.ascontiguousarray(targets, dtype="<f8")
    psi = np.ascontiguousarray(psi, dtype="<f8")
    _validate_arrays(targets, psi)
    batch_dir.mkdir(parents=True, exist_ok=True)
    _save_npy(batch_dir / "targets.npy", targets)
    _save_npy(batch_dir / "psi.npy", psi)
    path = batch_dir / "batch.nml"
    temporary = path.with_name("batch.nml.tmp")
    f90nml.Namelist({"vdsolverf_batch": config}).write(temporary, force=True)
    temporary.replace(path)
    return path


def read_batch_results(batch_dir) -> DustCurrentResult:
    """Collect all completed chunks in target order; reject missing/partial data.

    sample_count and terminal_box_count retain the output's int64 dtype.
    Files without a .done marker are never accepted. A summary is optional
    (all chunks can finish just before rank zero writes it).
    """
    path, config, targets, psi = _read_inputs(batch_dir)
    output = path.parent / config["output_dir"]
    arrays = {"flux": np.empty((len(targets), len(psi)), dtype=np.float64),
              "flux_error": np.empty((len(targets), len(psi)), dtype=np.float64),
              "status": np.empty(len(targets), dtype=np.int32),
              "sample_count": np.empty(len(targets), dtype=np.int64),
              "terminal_box_count": np.empty(len(targets), dtype=np.int64)}
    for start in range(0, len(targets), config["chunk_size"]):
        prefix = output / f"chunk_{start:09d}"
        if not prefix.with_suffix(".done").is_file():
            raise ValueError(f"chunk {start} is incomplete (missing .done marker)")
        count = min(config["chunk_size"], len(targets) - start)
        for name, values in arrays.items():
            dtype = "<f8" if name in ("flux", "flux_error") else "<i4" if name == "status" else "<i8"
            shape = (count, len(psi)) if name in ("flux", "flux_error") else (count,)
            values[start:start+count] = _load_npy(Path(f"{prefix}_{name}.npy"), dtype, shape)
    if np.any((arrays["status"] < 0) | (arrays["status"] > 5)):
        raise ValueError("invalid batch status")
    if any(np.any(arrays[name] < 0) for name in ("sample_count", "terminal_box_count")):
        raise ValueError("negative batch diagnostic counts")
    metadata = dict(config, batch_namelist=str(path), units="EMSES normalised; source density",
                    integration_rule="terminal_box_tensor_trapezoid")
    summary = output / "summary.json"
    if summary.is_file():
        metadata["summary"] = json.loads(summary.read_text(encoding="utf-8"))
    return DustCurrentResult(targets[:, :3], psi, metadata=metadata, **arrays)


def write_batch_results(batch_dir, result: DustCurrentResult):
    """Write results for a batch prepared by write_batch, including .done markers.

    Positions and psi must match its inputs exactly. This also permits fixtures
    with int64 diagnostic counts without running the native executable.
    """
    path, config, targets, psi = _read_inputs(batch_dir)
    if not isinstance(result, DustCurrentResult):
        raise TypeError("result must be a DustCurrentResult")
    if not np.array_equal(result.positions, targets[:, :3]) or not np.array_equal(result.psi_grid, psi):
        raise ValueError("result positions/psi do not match the prepared batch")
    if np.any((result.status < 0) | (result.status > 5)):
        raise ValueError("invalid batch status")
    if np.any(result.sample_count < 0) or np.any(result.terminal_box_count < 0):
        raise ValueError("negative batch diagnostic counts")
    output = path.parent / config["output_dir"]
    output.mkdir(parents=True, exist_ok=True)
    if (output / "summary.json").exists():
        (output / "summary.json").unlink()
    for start in range(0, len(targets), config["chunk_size"]):
        prefix = output / f"chunk_{start:09d}"
        marker = prefix.with_suffix(".done")
        if marker.exists():
            marker.unlink()
        for name, dtype in (("flux", "<f8"), ("flux_error", "<f8"), ("status", "<i4"),
                            ("sample_count", "<i8"), ("terminal_box_count", "<i8")):
            values = getattr(result, name)[start:start+config["chunk_size"]]
            _save_npy(Path(f"{prefix}_{name}.npy"), np.ascontiguousarray(values, dtype=dtype))
        temporary = marker.with_suffix(".done.tmp")
        temporary.touch()
        temporary.replace(marker)
    summary = dict(schema="vdsolverf.batch_summary", n_targets=len(targets), npsi=len(psi),
                   status_counts={str(i): int(np.count_nonzero(result.status == i)) for i in range(6)})
    temporary = output / "summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output / "summary.json")
    return output


_RANGE_DEFAULTS = dict(mode="velocity_range", dt=0.25, max_step=1000,
                       use_adaptive_dt=False, coverage_sigma=4.0, safety_factor=1.25,
                       source_samples_per_cell=1, minimum_count=4, accumulator_cache_size=20000,
                       max_probability_types=100, output_dir="out")
_RANGE_FLOATS = ("vx_min", "vx_max", "vy_min", "vy_max", "vz_min", "vz_max",
                 "weight_sum", "confidence")


def _validate_range(config):
    unknown = set(config) - set(_RANGE_DEFAULTS) - {"ispec", "bundle_dir", "n_threads"}
    if unknown:
        raise ValueError(f"unknown range batch settings: {sorted(unknown)}")
    if config["mode"] != "velocity_range":
        raise ValueError("expected mode='velocity_range'")
    for name in ("ispec", "max_step", "source_samples_per_cell", "minimum_count",
                 "accumulator_cache_size", "max_probability_types", "n_threads"):
        if name not in config:
            continue
        value = config[name]
        minimum = 0 if name == "ispec" else 1
        if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
                or not minimum <= value < np.iinfo(np.int32).max):
            raise ValueError(f"invalid range batch {name}")
        config[name] = int(value)
    for name in ("dt", "coverage_sigma", "safety_factor"):
        value = float(config[name])
        if not np.isfinite(value) or value <= 0 or (name == "safety_factor" and value < 1):
            raise ValueError(f"invalid range batch {name}")
        config[name] = value
    if not isinstance(config["use_adaptive_dt"], (bool, np.bool_)):
        raise ValueError("use_adaptive_dt must be a boolean")
    config["use_adaptive_dt"] = bool(config["use_adaptive_dt"])
    for name in ("bundle_dir", "output_dir"):
        value = os.fspath(config[name])
        if not value or len(value.encode()) > 4095 or any(c in value for c in "\x00\n\r"):
            raise ValueError(f"invalid range batch {name}")
        config[name] = value


def write_range_batch(batch_dir, bundle_dir, ispec, **settings):
    """Prepare a Python-free tier-2 run; return batch.nml (no targets/psi needed).

    Defaults match estimate_velocity_range_map. coverage_mode/eps_rel are
    resolved to coverage_sigma here. Bundle field flags are inherited; explicit
    use_*_field settings must agree. Omit n_threads for OMP_NUM_THREADS at runtime.
    """
    from ..emses.wrapper import _resolve_coverage_sigma

    batch_dir = Path(batch_dir).resolve()
    bundle = bundle_dir if isinstance(bundle_dir, FieldBundle) else read_field_bundle(bundle_dir)
    for name in ("use_electric_field", "use_magnetic_field"):
        if name in settings and settings.pop(name) != bundle.metadata[name]:
            raise ValueError(f"{name} must match the field bundle")
    config = dict(_RANGE_DEFAULTS)
    config.update(settings)
    config["coverage_sigma"] = _resolve_coverage_sigma(
        config["coverage_sigma"], coverage_mode=config.pop("coverage_mode", "sigma"),
        eps_rel=config.pop("eps_rel", 1e-6))
    if config.get("n_threads", 1) is None:
        config.pop("n_threads")
    config.update(ispec=ispec, bundle_dir=os.path.relpath(bundle.namelist_path.parent, batch_dir))
    _validate_range(config)
    if config["ispec"] not in bundle.metadata["ispecs"]:
        raise ValueError(f"species {ispec} is not present in the field bundle")
    batch_dir.mkdir(parents=True, exist_ok=True)
    path = batch_dir / "batch.nml"
    temporary = path.with_name("batch.nml.tmp")
    f90nml.Namelist({"vdsolverf_batch": config}).write(temporary, force=True)
    temporary.replace(path)
    return path


def read_range_batch_results(batch_dir) -> VelocityRangeMap:
    """Read a completed velocity_range job, with metadata and save() support.

    summary.json is the atomic completion marker. Moments are not collected.
    """
    path = _batch_path(batch_dir)
    config = dict(_RANGE_DEFAULTS)
    config.update(f90nml.read(path)["vdsolverf_batch"])
    _validate_range(config)
    output = path.parent / config["output_dir"] / "velocity_range"
    summary_path = output / "summary.json"
    if not summary_path.is_file():
        raise ValueError("velocity range batch is incomplete (missing summary.json)")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("schema") != "vdsolverf.velocity_range_summary" or summary.get("version") != "0.1.0":
        raise ValueError("unsupported velocity range summary schema/version")
    shape = tuple(summary["shape"])
    if len(shape) != 3 or any(type(n) is not int or n < 1 for n in shape):
        raise ValueError("invalid velocity range shape")
    arrays = {name: _load_npy(output / f"{name}.npy", "<f8", shape) for name in _RANGE_FLOATS}
    arrays.update({name: _load_npy(output / f"{name}.npy", "<i4", shape) for name in ("count", "status")})
    if np.any(arrays["count"] < 0) or np.any((arrays["status"] < 0) | (arrays["status"] > 3)):
        raise ValueError("invalid velocity range count/status")
    metadata = dict(config, **{k: summary[k] for k in ("directory", "istep")},
                    batch_namelist=str(path), field_bundle=str((path.parent / config["bundle_dir"]).resolve()),
                    collect_moments=False, show_progress=False, velocity_sample_mode="ellipsoid", summary=summary)
    metadata.update(summary["settings"])
    return VelocityRangeMap(x_edges=np.arange(shape[2]+1, dtype=np.float64),
                            y_edges=np.arange(shape[1]+1, dtype=np.float64),
                            z_edges=np.arange(shape[0]+1, dtype=np.float64),
                            metadata=metadata, directory=summary["directory"], **arrays)
