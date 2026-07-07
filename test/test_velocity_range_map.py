import unittest

import numpy as np

from vdsolverf.core import VelocityRangeMap


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


if __name__ == "__main__":
    unittest.main()
