import unittest

import numpy as np

from vdsolverf.core import Particle


class ParticleConstructorTest(unittest.TestCase):
    def test_accepts_documented_position_velocity_keywords(self):
        particle = Particle(position=[1.0, 2.0, 3.0], velocity=[4.0, 5.0, 6.0])

        np.testing.assert_allclose(particle.pos, [1.0, 2.0, 3.0])
        np.testing.assert_allclose(particle.vel, [4.0, 5.0, 6.0])

    def test_pos_vel_positional_constructor_still_works(self):
        particle = Particle([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])

        np.testing.assert_allclose(particle.pos, [1.0, 2.0, 3.0])
        np.testing.assert_allclose(particle.vel, [4.0, 5.0, 6.0])


if __name__ == "__main__":
    unittest.main()
