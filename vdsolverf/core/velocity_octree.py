from dataclasses import dataclass, field
from typing import Dict

import numpy as np


@dataclass
class VelocityOctreeResult:
    spatial_points: np.ndarray
    velocities: np.ndarray
    probabilities: np.ndarray
    spatial_index: np.ndarray
    leaf_spatial_index: np.ndarray
    leaf_bounds: np.ndarray
    leaf_value_min: np.ndarray
    leaf_value_max: np.ndarray
    leaf_depth: np.ndarray
    leaf_sample_start: np.ndarray
    leaf_sample_count: np.ndarray
    status: np.ndarray
    sample_count: np.ndarray
    leaf_count: np.ndarray
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        self.spatial_points = np.asarray(self.spatial_points, dtype=np.float64)
        self.velocities = np.asarray(self.velocities, dtype=np.float64)
        self.probabilities = np.asarray(self.probabilities, dtype=np.float64)
        self.spatial_index = np.asarray(self.spatial_index, dtype=np.int32)
        self.leaf_spatial_index = np.asarray(self.leaf_spatial_index, dtype=np.int32)
        self.leaf_bounds = np.asarray(self.leaf_bounds, dtype=np.float64)
        self.leaf_value_min = np.asarray(self.leaf_value_min, dtype=np.float64)
        self.leaf_value_max = np.asarray(self.leaf_value_max, dtype=np.float64)
        self.leaf_depth = np.asarray(self.leaf_depth, dtype=np.int32)
        self.leaf_sample_start = np.asarray(self.leaf_sample_start, dtype=np.int32)
        self.leaf_sample_count = np.asarray(self.leaf_sample_count, dtype=np.int32)
        self.status = np.asarray(self.status, dtype=np.int32)
        self.sample_count = np.asarray(self.sample_count, dtype=np.int32)
        self.leaf_count = np.asarray(self.leaf_count, dtype=np.int32)
        self.metadata = dict(self.metadata)

        self._validate_shapes()

    @property
    def nspatial(self) -> int:
        return int(self.spatial_points.shape[0])

    @property
    def nsamples(self) -> int:
        return int(self.probabilities.shape[0])

    @property
    def nleaves(self) -> int:
        return int(self.leaf_bounds.shape[0])

    @property
    def valid_probability_mask(self) -> np.ndarray:
        return np.isfinite(self.probabilities)

    def _validate_shapes(self):
        if self.spatial_points.ndim != 2 or self.spatial_points.shape[1] != 3:
            raise ValueError("spatial_points must have shape (nspatial, 3)")
        if self.velocities.ndim != 2 or self.velocities.shape[1] != 3:
            raise ValueError("velocities must have shape (nsample, 3)")
        if self.probabilities.shape != (self.velocities.shape[0],):
            raise ValueError("probabilities must have shape (nsample,)")
        if self.spatial_index.shape != (self.velocities.shape[0],):
            raise ValueError("spatial_index must have shape (nsample,)")
        if self.leaf_bounds.ndim != 2 or self.leaf_bounds.shape[1] != 6:
            raise ValueError("leaf_bounds must have shape (nleaf, 6)")

        nleaf = self.leaf_bounds.shape[0]
        for name in (
            "leaf_spatial_index",
            "leaf_value_min",
            "leaf_value_max",
            "leaf_depth",
            "leaf_sample_start",
            "leaf_sample_count",
        ):
            if getattr(self, name).shape != (nleaf,):
                raise ValueError(f"{name} must have shape (nleaf,)")

        nspatial = self.spatial_points.shape[0]
        for name in ("status", "sample_count", "leaf_count"):
            if getattr(self, name).shape != (nspatial,):
                raise ValueError(f"{name} must have shape (nspatial,)")
