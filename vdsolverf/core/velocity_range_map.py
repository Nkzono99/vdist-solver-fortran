from dataclasses import dataclass, field
from itertools import product
from typing import Dict, List, Tuple

import numpy as np

from .particles import Particle


VelocityBins = Tuple[int, int, int]


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

        self._validate_shapes()

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
