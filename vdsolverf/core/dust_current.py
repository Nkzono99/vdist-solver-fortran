"""Dust OML flux moments; multiply flux by species charge and density."""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class DustCurrentResult:
    positions: np.ndarray
    psi_grid: np.ndarray
    flux: np.ndarray
    flux_error: np.ndarray
    status: np.ndarray
    sample_count: np.ndarray
    terminal_box_count: np.ndarray
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in ("positions", "psi_grid", "flux", "flux_error"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=np.float64))
        for name in ("status", "sample_count", "terminal_box_count"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=np.int32))
        self.metadata = dict(self.metadata)
        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError("positions must have shape (nspatial, 3)")
        if self.psi_grid.ndim != 1:
            raise ValueError("psi_grid must have shape (npsi,)")
        expected = (len(self.positions), len(self.psi_grid))
        if self.flux.shape != expected or self.flux_error.shape != expected:
            raise ValueError("flux and flux_error must have shape (nspatial, npsi)")
        for name in ("status", "sample_count", "terminal_box_count"):
            if getattr(self, name).shape != (len(self.positions),):
                raise ValueError(f"{name} must have shape (nspatial,)")
