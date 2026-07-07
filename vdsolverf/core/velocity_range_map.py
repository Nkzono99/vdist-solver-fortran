from dataclasses import dataclass, field
from itertools import product
import json
from pathlib import Path
from typing import Dict, List, Tuple, Union

import numpy as np

from .particles import Particle


VelocityBins = Tuple[int, int, int]
PathLike = Union[str, Path]


@dataclass
class VelocityRangeIndex:
    shape: Tuple[int, int, int, int, int, int]
    flat_indices: np.ndarray

    def reshape(self, values: np.ndarray, fill_value: float = np.nan) -> np.ndarray:
        values = np.asarray(values)
        if values.shape[0] != self.flat_indices.shape[0]:
            raise ValueError(
                "values length must match the number of particles created by the index"
            )

        result_dtype = np.result_type(values.dtype, type(fill_value))
        result = np.full(self.shape, fill_value, dtype=result_dtype)
        result.reshape(-1)[self.flat_indices] = values
        return result


@dataclass(frozen=True)
class VelocityRangeCell:
    range_map: "VelocityRangeMap"
    iz: int
    iy: int
    ix: int

    @property
    def index(self) -> Tuple[int, int, int]:
        return self.iz, self.iy, self.ix

    @property
    def position(self) -> np.ndarray:
        return np.array(
            [
                0.5 * (self.range_map.x_edges[self.ix] + self.range_map.x_edges[self.ix + 1]),
                0.5 * (self.range_map.y_edges[self.iy] + self.range_map.y_edges[self.iy + 1]),
                0.5 * (self.range_map.z_edges[self.iz] + self.range_map.z_edges[self.iz + 1]),
            ],
            dtype=np.float64,
        )

    @property
    def vmin(self) -> np.ndarray:
        return np.array(
            [
                self.range_map.vx_min[self.index],
                self.range_map.vy_min[self.index],
                self.range_map.vz_min[self.index],
            ],
            dtype=np.float64,
        )

    @property
    def vmax(self) -> np.ndarray:
        return np.array(
            [
                self.range_map.vx_max[self.index],
                self.range_map.vy_max[self.index],
                self.range_map.vz_max[self.index],
            ],
            dtype=np.float64,
        )

    @property
    def count(self) -> int:
        return int(self.range_map.count[self.index])

    @property
    def weight_sum(self) -> float:
        return float(self.range_map.weight_sum[self.index])

    @property
    def mean_v(self) -> np.ndarray:
        return self.range_map.mean_v[self.index]

    @property
    def cov_v(self) -> np.ndarray:
        return self.range_map.cov_v[self.index]

    @property
    def status(self) -> int:
        return int(self.range_map.status[self.index])

    @property
    def confidence(self) -> float:
        return float(self.range_map.confidence[self.index])

    @property
    def valid(self) -> bool:
        return bool(self.range_map.valid_mask[self.index])

    def velocity_axes(self, velocity_bins: VelocityBins) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        nvx, nvy, nvz = velocity_bins
        if min(nvx, nvy, nvz) < 1:
            raise ValueError("velocity_bins must contain positive integers")

        return (
            np.linspace(self.range_map.vx_min[self.index], self.range_map.vx_max[self.index], nvx),
            np.linspace(self.range_map.vy_min[self.index], self.range_map.vy_max[self.index], nvy),
            np.linspace(self.range_map.vz_min[self.index], self.range_map.vz_max[self.index], nvz),
        )

    def create_particles(
        self,
        velocity_bins: VelocityBins,
        *,
        include_invalid: bool = False,
    ) -> Tuple[List[Particle], VelocityRangeIndex]:
        nvx, nvy, nvz = velocity_bins
        if min(nvx, nvy, nvz) < 1:
            raise ValueError("velocity_bins must contain positive integers")

        shape = (nvz, nvy, nvx)
        if not include_invalid and not self.valid:
            return [], VelocityRangeIndex(shape, np.array([], dtype=np.int64))

        vx_values, vy_values, vz_values = self.velocity_axes(velocity_bins)
        position = self.position
        particles: List[Particle] = []
        flat_indices = []

        for ivz, ivy, ivx in product(range(nvz), range(nvy), range(nvx)):
            particles.append(
                Particle(
                    position.copy(),
                    np.array(
                        [vx_values[ivx], vy_values[ivy], vz_values[ivz]],
                        dtype=np.float64,
                    ),
                )
            )
            flat_indices.append(np.ravel_multi_index((ivz, ivy, ivx), shape))

        return particles, VelocityRangeIndex(shape, np.array(flat_indices, dtype=np.int64))


@dataclass
class VelocityRangeMap:
    x_edges: np.ndarray
    y_edges: np.ndarray
    z_edges: np.ndarray

    vx_min: np.ndarray
    vx_max: np.ndarray
    vy_min: np.ndarray
    vy_max: np.ndarray
    vz_min: np.ndarray
    vz_max: np.ndarray

    count: np.ndarray
    weight_sum: np.ndarray = None
    mean_v: np.ndarray = None
    cov_v: np.ndarray = None
    status: np.ndarray = None
    confidence: np.ndarray = None
    metadata: Dict = field(default_factory=dict)
    directory: PathLike = None

    def __post_init__(self):
        self.x_edges = np.asarray(self.x_edges, dtype=np.float64)
        self.y_edges = np.asarray(self.y_edges, dtype=np.float64)
        self.z_edges = np.asarray(self.z_edges, dtype=np.float64)

        for name in (
            "vx_min",
            "vx_max",
            "vy_min",
            "vy_max",
            "vz_min",
            "vz_max",
        ):
            setattr(self, name, np.asarray(getattr(self, name), dtype=np.float64))

        self.count = np.asarray(self.count, dtype=np.int32)
        cell_shape = self.count.shape

        if self.weight_sum is None:
            self.weight_sum = np.zeros(cell_shape, dtype=np.float64)
        else:
            self.weight_sum = np.asarray(self.weight_sum, dtype=np.float64)

        if self.mean_v is None:
            self.mean_v = np.full(cell_shape + (3,), np.nan, dtype=np.float64)
        else:
            self.mean_v = np.asarray(self.mean_v, dtype=np.float64)

        if self.cov_v is None:
            self.cov_v = np.full(cell_shape + (3, 3), np.nan, dtype=np.float64)
        else:
            self.cov_v = np.asarray(self.cov_v, dtype=np.float64)

        if self.status is None:
            self.status = np.where(self.count > 0, 0, 2).astype(np.int32)
        else:
            self.status = np.asarray(self.status, dtype=np.int32)

        if self.confidence is None:
            self.confidence = np.zeros(cell_shape, dtype=np.float64)
        else:
            self.confidence = np.asarray(self.confidence, dtype=np.float64)

        self.metadata = dict(self.metadata)
        if self.directory is not None:
            self.directory = Path(self.directory)

        self._validate_shapes()

    def __getitem__(self, key) -> VelocityRangeCell:
        if not isinstance(key, tuple) or len(key) != 3:
            raise IndexError("VelocityRangeMap expects indices as range_map[iz, iy, ix]")

        iz, iy, ix = (
            self._normalize_index(key[0], self.cell_shape[0], "z"),
            self._normalize_index(key[1], self.cell_shape[1], "y"),
            self._normalize_index(key[2], self.cell_shape[2], "x"),
        )
        return VelocityRangeCell(self, iz, iy, ix)

    @property
    def cell_shape(self) -> Tuple[int, int, int]:
        return self.count.shape

    @property
    def valid_mask(self) -> np.ndarray:
        return (
            (self.count > 0)
            & np.isfinite(self.vx_min)
            & np.isfinite(self.vx_max)
            & np.isfinite(self.vy_min)
            & np.isfinite(self.vy_max)
            & np.isfinite(self.vz_min)
            & np.isfinite(self.vz_max)
            & (self.vx_max >= self.vx_min)
            & (self.vy_max >= self.vy_min)
            & (self.vz_max >= self.vz_min)
        )

    def create_particles(
        self,
        velocity_bins: VelocityBins,
        *,
        include_invalid: bool = False,
    ) -> Tuple[List[Particle], VelocityRangeIndex]:
        nvx, nvy, nvz = velocity_bins
        if min(nvx, nvy, nvz) < 1:
            raise ValueError("velocity_bins must contain positive integers")

        particles: List[Particle] = []
        flat_indices = []
        nz, ny, nx = self.cell_shape
        shape = (nz, ny, nx, nvz, nvy, nvx)
        valid = self.valid_mask

        x_centers = 0.5 * (self.x_edges[:-1] + self.x_edges[1:])
        y_centers = 0.5 * (self.y_edges[:-1] + self.y_edges[1:])
        z_centers = 0.5 * (self.z_edges[:-1] + self.z_edges[1:])

        for iz, iy, ix in product(range(nz), range(ny), range(nx)):
            if not include_invalid and not valid[iz, iy, ix]:
                continue

            vx_values = np.linspace(self.vx_min[iz, iy, ix], self.vx_max[iz, iy, ix], nvx)
            vy_values = np.linspace(self.vy_min[iz, iy, ix], self.vy_max[iz, iy, ix], nvy)
            vz_values = np.linspace(self.vz_min[iz, iy, ix], self.vz_max[iz, iy, ix], nvz)
            position = np.array(
                [x_centers[ix], y_centers[iy], z_centers[iz]],
                dtype=np.float64,
            )

            for ivz, ivy, ivx in product(range(nvz), range(nvy), range(nvx)):
                particles.append(
                    Particle(
                        position.copy(),
                        np.array(
                            [vx_values[ivx], vy_values[ivy], vz_values[ivz]],
                            dtype=np.float64,
                        ),
                    )
                )
                flat_indices.append(np.ravel_multi_index((iz, iy, ix, ivz, ivy, ivx), shape))

        return particles, VelocityRangeIndex(shape, np.array(flat_indices, dtype=np.int64))

    def default_path(self, filename: str = None) -> Path:
        if self.directory is None:
            raise ValueError("directory is required to build the default range-map path")

        return self.default_path_for(
            self.directory,
            ispec=self.metadata.get("ispec"),
            istep=self.metadata.get("istep"),
            filename=filename,
        )

    def save(
        self,
        path: PathLike = None,
        *,
        directory: PathLike = None,
        filename: str = None,
    ) -> Path:
        if path is not None and (directory is not None or filename is not None):
            raise ValueError("pass either path or directory/filename, not both")

        stored_directory = self.directory
        if path is None:
            save_directory = Path(directory) if directory is not None else self.directory
            if save_directory is None:
                raise ValueError("path or directory is required to save a range map")
            path = self.default_path_for(
                save_directory,
                ispec=self.metadata.get("ispec"),
                istep=self.metadata.get("istep"),
                filename=filename,
            )
            if stored_directory is None:
                stored_directory = save_directory
        else:
            path = Path(path)

        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            x_edges=self.x_edges,
            y_edges=self.y_edges,
            z_edges=self.z_edges,
            vx_min=self.vx_min,
            vx_max=self.vx_max,
            vy_min=self.vy_min,
            vy_max=self.vy_max,
            vz_min=self.vz_min,
            vz_max=self.vz_max,
            count=self.count,
            weight_sum=self.weight_sum,
            mean_v=self.mean_v,
            cov_v=self.cov_v,
            status=self.status,
            confidence=self.confidence,
            metadata_json=np.array(json.dumps(_json_ready(self.metadata))),
            directory=np.array("" if stored_directory is None else str(stored_directory)),
        )
        return path

    @classmethod
    def load(
        cls,
        path: PathLike = None,
        *,
        directory: PathLike = None,
        ispec: int = None,
        istep: int = None,
        filename: str = None,
    ) -> "VelocityRangeMap":
        if path is not None and (directory is not None or filename is not None):
            raise ValueError("pass either path or directory/filename, not both")

        load_directory = None if directory is None else Path(directory)
        if path is None:
            if load_directory is None:
                raise ValueError("path or directory is required to load a range map")
            path = cls.default_path_for(
                load_directory,
                ispec=ispec,
                istep=istep,
                filename=filename,
            )
        else:
            path = Path(path)

        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"].item()))
            stored_directory = str(data["directory"].item())
            map_directory = Path(stored_directory) if stored_directory else load_directory

            return cls(
                x_edges=data["x_edges"],
                y_edges=data["y_edges"],
                z_edges=data["z_edges"],
                vx_min=data["vx_min"],
                vx_max=data["vx_max"],
                vy_min=data["vy_min"],
                vy_max=data["vy_max"],
                vz_min=data["vz_min"],
                vz_max=data["vz_max"],
                count=data["count"],
                weight_sum=data["weight_sum"],
                mean_v=data["mean_v"],
                cov_v=data["cov_v"],
                status=data["status"],
                confidence=data["confidence"],
                metadata=metadata,
                directory=map_directory,
            )

    @staticmethod
    def default_filename_for(ispec: int = None, istep: int = None) -> str:
        if ispec is not None and istep is not None:
            return f"vdsolverf-velocity-range-map-ispec{int(ispec)}-istep{int(istep)}.npz"
        return "vdsolverf-velocity-range-map.npz"

    @classmethod
    def default_path_for(
        cls,
        directory: PathLike,
        *,
        ispec: int = None,
        istep: int = None,
        filename: str = None,
    ) -> Path:
        return Path(directory) / (filename or cls.default_filename_for(ispec, istep))

    def expand_cells(self, mask: np.ndarray, factor: float):
        if factor <= 0:
            raise ValueError("factor must be positive")

        mask = np.asarray(mask, dtype=bool)
        if mask.shape != self.cell_shape:
            raise ValueError("mask shape must match range map cells")

        for vmin_name, vmax_name in (
            ("vx_min", "vx_max"),
            ("vy_min", "vy_max"),
            ("vz_min", "vz_max"),
        ):
            vmin = getattr(self, vmin_name)
            vmax = getattr(self, vmax_name)
            center = 0.5 * (vmin + vmax)
            half_width = 0.5 * (vmax - vmin) * factor
            vmin[mask] = center[mask] - half_width[mask]
            vmax[mask] = center[mask] + half_width[mask]

    def _validate_shapes(self):
        expected = (
            len(self.z_edges) - 1,
            len(self.y_edges) - 1,
            len(self.x_edges) - 1,
        )
        if self.count.shape != expected:
            raise ValueError("count shape must match z/y/x edge dimensions")

        for name in (
            "vx_min",
            "vx_max",
            "vy_min",
            "vy_max",
            "vz_min",
            "vz_max",
            "weight_sum",
            "status",
            "confidence",
        ):
            if getattr(self, name).shape != expected:
                raise ValueError(f"{name} shape must match count shape")

        if self.mean_v.shape != expected + (3,):
            raise ValueError("mean_v shape must be count.shape + (3,)")
        if self.cov_v.shape != expected + (3, 3):
            raise ValueError("cov_v shape must be count.shape + (3, 3)")

    @staticmethod
    def _normalize_index(value, size: int, axis_name: str) -> int:
        if isinstance(value, slice):
            raise IndexError("VelocityRangeMap cell access does not support slices")

        index = int(value)
        if index < 0:
            index += size
        if index < 0 or index >= size:
            raise IndexError(f"{axis_name} index out of range")
        return index


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value
