import unittest

import numpy as np

from vdsolverf.core import VelocityOctreeResult


class VelocityOctreeResultTest(unittest.TestCase):
    def test_stores_ragged_samples_and_leaf_diagnostics(self):
        result = VelocityOctreeResult(
            spatial_points=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
            velocities=np.array([[10.0, 0.0, 0.0], [0.0, 20.0, 0.0]]),
            probabilities=np.array([0.5, np.nan]),
            spatial_index=np.array([0, 1], dtype=np.int32),
            leaf_spatial_index=np.array([0], dtype=np.int32),
            leaf_bounds=np.array([[-1.0, 1.0, -2.0, 2.0, -3.0, 3.0]]),
            leaf_value_min=np.array([0.1]),
            leaf_value_max=np.array([0.5]),
            leaf_depth=np.array([2], dtype=np.int32),
            leaf_sample_start=np.array([0], dtype=np.int32),
            leaf_sample_count=np.array([2], dtype=np.int32),
            status=np.array([0, 1], dtype=np.int32),
            sample_count=np.array([1, 1], dtype=np.int32),
            leaf_count=np.array([1, 0], dtype=np.int32),
            metadata={"max_depth": 3},
        )

        self.assertEqual(result.nspatial, 2)
        self.assertEqual(result.nsamples, 2)
        self.assertEqual(result.nleaves, 1)
        np.testing.assert_array_equal(result.valid_probability_mask, [True, False])
        self.assertEqual(result.metadata["max_depth"], 3)


if __name__ == "__main__":
    unittest.main()
