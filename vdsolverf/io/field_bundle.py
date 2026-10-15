"""Portable, reusable relocated EMSES fields for ctypes and standalone batches."""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib.metadata import version as distribution_version
import json
from pathlib import Path
import shutil
import uuid

import numpy as np


@dataclass
class FieldBundle:
    namelist_path: Path
    metadata: dict
    _arrays: dict = field(default_factory=dict, init=False, repr=False)

    def get_ebvalues(self, ispec: int) -> np.ndarray:
        """Lazily memory-map a species' little-endian, C-order float64 fields."""
        if ispec not in self.metadata["ispecs"]:
            raise ValueError(f"species {ispec} is not present in the field bundle")
        if ispec not in self._arrays:
            path = self.namelist_path.parent / f"ebvalues_ispec{ispec}.npy"
            with path.open("rb") as stream:
                if np.lib.format.read_magic(stream) != (1, 0):
                    raise ValueError("field bundle requires npy v1.0")
            values = np.load(path, mmap_mode="r", allow_pickle=False)
            expected = tuple(self.metadata[k] + 1 for k in ("nz", "ny", "nx")) + (9,)
            if values.shape != expected or values.dtype.str != "<f8" or not values.flags.c_contiguous:
                raise ValueError(f"invalid field array {path}: expected C-order <f8 {expected}")
            self._arrays[ispec] = values
        return self._arrays[ispec]

    def __getitem__(self, ispec):
        return self.get_ebvalues(ispec)


def read_field_bundle(path) -> FieldBundle:
    """Read bundle.json from a bundle directory or a JSON path; fields stay lazy."""
    path = Path(path).resolve()
    if path.is_dir():
        path = path / "bundle.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "vdsolverf.field_bundle" or metadata.get("version") != "0.1.0":
        raise ValueError("unsupported field bundle schema/version")
    for key in ("nx", "ny", "nz"):
        if type(metadata.get(key)) is not int or metadata[key] < 1:
            raise ValueError(f"invalid field bundle {key}")
    species = metadata.get("ispecs")
    if not isinstance(species, list) or not species or any(type(s) is not int or s < 0 for s in species):
        raise ValueError("invalid field bundle ispecs")
    if len(set(species)) != len(species):
        raise ValueError("duplicate field bundle species")
    for key in ("use_electric_field", "use_magnetic_field"):
        if type(metadata.get(key)) is not bool:
            raise ValueError(f"invalid field bundle {key}")
    if type(metadata.get("istep")) is not int or not isinstance(metadata.get("directory"), str):
        raise ValueError("invalid field bundle directory/istep")
    namelist = path.parent / "plasma-vdsolverf.inp"
    if not namelist.is_file():
        raise ValueError(f"field bundle namelist is missing: {namelist}")
    return FieldBundle(namelist, metadata)


def write_field_bundle(directory, istep, ispecs, output_dir, *,
                       use_electric_field=True, use_magnetic_field=True) -> FieldBundle:
    """Write exactly TempolaryInput's namelist and create_relocated_ebvalues arrays.

    Species are zero-based, matching the Python EMSES APIs. Npy files explicitly
    use version 1.0, dtype <f8 and C order for the standalone reader.
    """
    # Import lazily to keep io usable without importing ctypes or loading fields.
    from ..emses import wrapper

    species = list(ispecs)
    if not species or any(not isinstance(s, (int, np.integer)) or s < 0 for s in species):
        raise ValueError("ispecs must contain non-negative integer species indices")
    species = [int(s) for s in species]
    if len(set(species)) != len(species):
        raise ValueError("ispecs must not contain duplicates")
    data = wrapper.emout.Emout(directory)
    if any(s >= data.inp.nspec for s in species):
        raise ValueError("ispecs contains an index outside the namelist species range")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    # A unique source filename also permits output_dir == directory safely.
    with wrapper.TempolaryInput(data, suffix=f"bundle-{uuid.uuid4().hex}") as tmp:
        shutil.copyfile(tmp.tmppath, output_dir / "plasma-vdsolverf.inp")
    for ispec in species:
        values = wrapper.create_relocated_ebvalues(
            data, istep, ispec=ispec, use_electric_field=use_electric_field,
            use_magnetic_field=use_magnetic_field,
        )
        values = np.ascontiguousarray(values, dtype="<f8")
        with (output_dir / f"ebvalues_ispec{ispec}.npy").open("wb") as stream:
            np.lib.format.write_array(stream, values, version=(1, 0), allow_pickle=False)
    metadata = {
        "schema": "vdsolverf.field_bundle", "version": "0.1.0",
        "directory": str(Path(directory).resolve()), "istep": int(istep),
        "nx": int(data.inp.nx), "ny": int(data.inp.ny), "nz": int(data.inp.nz),
        "ispecs": species, "use_electric_field": bool(use_electric_field),
        "use_magnetic_field": bool(use_magnetic_field),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "vdsolverf_version": distribution_version("vdist-solver-fortran"),
    }
    (output_dir / "bundle.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return read_field_bundle(output_dir)
