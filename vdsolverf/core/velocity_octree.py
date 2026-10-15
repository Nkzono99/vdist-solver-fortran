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
    leaf_is_terminal: np.ndarray = None

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
        if self.leaf_is_terminal is None:
            # Compatibility with saved/constructed results and older libraries:
            # dyadic box addresses identify parents without a quadratic search.
            self.leaf_is_terminal = np.ones(self.nleaves, dtype=bool)
            for spatial in np.unique(self.leaf_spatial_index):
                indexes = np.flatnonzero(self.leaf_spatial_index == spatial)
                roots = indexes[self.leaf_depth[indexes] == 0]
                if not roots.size:
                    continue
                root = self.leaf_bounds[roots[0]]
                widths = root[1::2] - root[::2]
                nodes = {}
                for i in indexes:
                    depth = int(self.leaf_depth[i])
                    center = (self.leaf_bounds[i, ::2] + self.leaf_bounds[i, 1::2]) * 0.5
                    address = tuple(np.floor((center - root[::2]) / widths * 2**depth).astype(int))
                    nodes[depth, address] = i
                for (depth, address) in nodes:
                    parent = nodes.get((depth - 1, tuple(a // 2 for a in address)))
                    if parent is not None:
                        self.leaf_is_terminal[parent] = False
        else:
            self.leaf_is_terminal = np.asarray(self.leaf_is_terminal, dtype=bool)
            if self.leaf_is_terminal.shape != (self.nleaves,):
                raise ValueError("leaf_is_terminal must have shape (nleaf,)")

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

    def sample_mask(self, spatial_index: int) -> np.ndarray:
        self._validate_spatial_index(spatial_index)
        return self.spatial_index == spatial_index

    def leaf_mask(self, spatial_index: int) -> np.ndarray:
        self._validate_spatial_index(spatial_index)
        return self.leaf_spatial_index == spatial_index

    def velocities_for_spatial(self, spatial_index: int) -> np.ndarray:
        return self.velocities[self.sample_mask(spatial_index)]

    def probabilities_for_spatial(self, spatial_index: int) -> np.ndarray:
        return self.probabilities[self.sample_mask(spatial_index)]

    def leaf_sample_slice(self, leaf_index: int) -> slice:
        if leaf_index < 0 or leaf_index >= self.nleaves:
            raise IndexError("leaf_index out of range")
        start = int(self.leaf_sample_start[leaf_index])
        stop = start + int(self.leaf_sample_count[leaf_index])
        return slice(start, stop)

    def _validate_spatial_index(self, spatial_index: int):
        if spatial_index < 0 or spatial_index >= self.nspatial:
            raise IndexError("spatial_index out of range")

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

        if np.any(self.spatial_index < 0) or np.any(self.spatial_index >= nspatial):
            raise ValueError("spatial_index contains out-of-range values")
        if np.any(self.leaf_spatial_index < 0) or np.any(self.leaf_spatial_index >= nspatial):
            raise ValueError("leaf_spatial_index contains out-of-range values")
        if np.any(self.sample_count < 0):
            raise ValueError("sample_count must be non-negative")
        if np.any(self.leaf_count < 0):
            raise ValueError("leaf_count must be non-negative")
        if int(np.sum(self.sample_count)) != self.nsamples:
            raise ValueError("sum(sample_count) must equal nsamples")
        if int(np.sum(self.leaf_count)) != self.nleaves:
            raise ValueError("sum(leaf_count) must equal nleaves")
        if np.any(self.leaf_sample_start < 0) or np.any(self.leaf_sample_count < 0):
            raise ValueError("leaf sample slices must be non-negative")
        if np.any(self.leaf_sample_start + self.leaf_sample_count > self.nsamples):
            raise ValueError("leaf sample slices must stay within compact samples")
