import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from vdsolverf.core import VelocityRangeCell, VelocityRangeMap


class VelocityRangeMapCreateParticlesTest(unittest.TestCase):
    def test_create_particles_returns_particles_and_reshape_index(self):
        range_map = VelocityRangeMap(
            x_edges=np.array([0.0, 1.0, 2.0]),
            y_edges=np.array([0.0, 1.0]),
            z_edges=np.array([0.0, 1.0]),
            vx_min=np.array([[[10.0, 0.0]]]),
            vx_max=np.array([[[12.0, 0.0]]]),
            vy_min=np.array([[[20.0, 0.0]]]),
            vy_max=np.array([[[20.0, 0.0]]]),
            vz_min=np.array([[[30.0, 0.0]]]),
            vz_max=np.array([[[32.0, 0.0]]]),
            count=np.array([[[4, 0]]], dtype=np.int32),
        )

        particles, index = range_map.create_particles(velocity_bins=(3, 1, 2))

        self.assertEqual(len(particles), 6)
        np.testing.assert_allclose(
            [particle.pos for particle in particles],
            [[0.5, 0.5, 0.5]] * 6,
        )
        np.testing.assert_allclose(
            [particle.vel for particle in particles],
            [
                [10.0, 20.0, 30.0],
                [11.0, 20.0, 30.0],
                [12.0, 20.0, 30.0],
                [10.0, 20.0, 32.0],
                [11.0, 20.0, 32.0],
                [12.0, 20.0, 32.0],
            ],
        )

        reshaped = index.reshape(np.arange(6, dtype=np.float64))

        self.assertEqual(reshaped.shape, (1, 1, 2, 2, 1, 3))
        np.testing.assert_allclose(reshaped[0, 0, 0, :, :, :].ravel(), np.arange(6))
        self.assertTrue(np.isnan(reshaped[0, 0, 1]).all())

    def test_cell_view_creates_particles_for_one_cell(self):
        range_map = sample_range_map()

        cell = range_map[0, 0, 1]

        self.assertIsInstance(cell, VelocityRangeCell)
        self.assertEqual(cell.index, (0, 0, 1))
        self.assertTrue(cell.valid)
        np.testing.assert_allclose(cell.position, [1.5, 0.5, 0.5])
        np.testing.assert_allclose(cell.vmin, [1.0, 10.0, 100.0])
        np.testing.assert_allclose(cell.vmax, [3.0, 12.0, 102.0])
        self.assertEqual(cell.count, 5)

        vx, vy, vz = cell.velocity_axes((3, 2, 2))
        np.testing.assert_allclose(vx, [1.0, 2.0, 3.0])
        np.testing.assert_allclose(vy, [10.0, 12.0])
        np.testing.assert_allclose(vz, [100.0, 102.0])

        particles, index = cell.create_particles(velocity_bins=(3, 2, 2))

        self.assertEqual(len(particles), 12)
        self.assertEqual(index.shape, (2, 2, 3))
        np.testing.assert_allclose(particles[0].pos, [1.5, 0.5, 0.5])
        np.testing.assert_allclose(particles[0].vel, [1.0, 10.0, 100.0])
        np.testing.assert_allclose(particles[-1].vel, [3.0, 12.0, 102.0])
        np.testing.assert_allclose(index.reshape(np.arange(12)).ravel(), np.arange(12))

    def test_invalid_cell_view_returns_no_particles_by_default(self):
        range_map = sample_range_map()

        particles, index = range_map[0, 0, 0].create_particles(velocity_bins=(2, 2, 2))

        self.assertEqual(particles, [])
        self.assertEqual(index.shape, (2, 2, 2))
        self.assertTrue(np.isnan(index.reshape(np.array([], dtype=np.float64))).all())

    def test_save_and_load_default_path(self):
        with TemporaryDirectory() as tmpdir:
            range_map = sample_range_map(
                directory=tmpdir,
                metadata={"ispec": 0, "istep": -1, "coverage_sigma": 4.0},
            )

            saved_path = range_map.save()

            self.assertEqual(
                saved_path,
                Path(tmpdir) / "vdsolverf-velocity-range-map-ispec0-istep-1.npz",
            )
            self.assertTrue(saved_path.exists())

            loaded = VelocityRangeMap.load(directory=tmpdir, ispec=0, istep=-1)

            self.assertEqual(loaded.directory, Path(tmpdir))
            self.assertEqual(loaded.metadata["ispec"], 0)
            self.assertEqual(loaded.metadata["istep"], -1)
            np.testing.assert_allclose(loaded.vx_min, range_map.vx_min)
            np.testing.assert_allclose(loaded.vx_max, range_map.vx_max)
            np.testing.assert_array_equal(loaded.count, range_map.count)

    def test_save_requires_path_or_directory(self):
        range_map = sample_range_map()

        with self.assertRaises(ValueError):
            range_map.save()


def sample_range_map(directory=None, metadata=None):
    return VelocityRangeMap(
        x_edges=np.array([0.0, 1.0, 2.0]),
        y_edges=np.array([0.0, 1.0]),
        z_edges=np.array([0.0, 1.0]),
        vx_min=np.array([[[np.nan, 1.0]]]),
        vx_max=np.array([[[np.nan, 3.0]]]),
        vy_min=np.array([[[np.nan, 10.0]]]),
        vy_max=np.array([[[np.nan, 12.0]]]),
        vz_min=np.array([[[np.nan, 100.0]]]),
        vz_max=np.array([[[np.nan, 102.0]]]),
        count=np.array([[[0, 5]]], dtype=np.int32),
        weight_sum=np.array([[[0.0, 5.0]]]),
        mean_v=np.array([[[[np.nan, np.nan, np.nan], [2.0, 11.0, 101.0]]]]),
        cov_v=np.array(
            [[[[[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]],
               [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]]]]
        ),
        status=np.array([[[2, 0]]], dtype=np.int32),
        confidence=np.array([[[0.0, 1.0]]]),
        metadata=metadata or {},
        directory=directory,
    )


if __name__ == "__main__":
    unittest.main()
