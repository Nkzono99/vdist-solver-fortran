import unittest

import numpy as np

from vdsolverf.core import PhaseGrid


class PhaseGridCreateParticlesTest(unittest.TestCase):
    def test_create_particles_does_not_materialize_full_grid(self):
        phase_grid = PhaseGrid(
            x=(1.0, 2.0, 2),
            y=3.0,
            z=4.0,
            vx=(-1.0, 1.0, 3),
            vy=5.0,
            vz=6.0,
        )

        def fail_if_full_grid_is_built():
            raise AssertionError("create_particles should not materialize create_grid()")

        phase_grid.create_grid = fail_if_full_grid_is_built

        particles = phase_grid.create_particles()

        self.assertEqual(len(particles), 6)
        np.testing.assert_allclose(
            [particle.pos for particle in particles],
            [
                [1.0, 3.0, 4.0],
                [1.0, 3.0, 4.0],
                [1.0, 3.0, 4.0],
                [2.0, 3.0, 4.0],
                [2.0, 3.0, 4.0],
                [2.0, 3.0, 4.0],
            ],
        )
        np.testing.assert_allclose(
            [particle.vel for particle in particles],
            [
                [-1.0, 5.0, 6.0],
                [0.0, 5.0, 6.0],
                [1.0, 5.0, 6.0],
                [-1.0, 5.0, 6.0],
                [0.0, 5.0, 6.0],
                [1.0, 5.0, 6.0],
            ],
        )


if __name__ == "__main__":
    unittest.main()
