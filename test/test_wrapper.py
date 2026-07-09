import inspect
import tempfile
import unittest
from pathlib import Path
from typing import Tuple, get_type_hints
from unittest.mock import patch

import numpy as np

from vdsolverf.core import Particle, PhaseGrid
from vdsolverf.emses import wrapper


class _FakeInput:
    nx = 1
    ny = 1
    nz = 1


class _FakeData:
    inp = _FakeInput()
    directory = Path("/tmp/vdsolverf-fake-run")


class _FakeTemporaryInput:
    def __init__(self, data, suffix=None):
        self.tmppath = Path(tempfile.gettempdir()) / "vdsolverf-test.inp"

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False


class _FakeBacktracesFunction:
    def __init__(self):
        self.argtypes = None
        self.restype = None
        self.n_threads = None

    def __call__(self, *args):
        self.n_threads = args[-1]._obj.value
        args[16][:] = 0.0
        args[17][:] = 0.0
        args[18][:] = 0.0
        args[19][:] = 0


class _FakeEstimateVelocityRangeMapFunction:
    def __init__(self):
        self.argtypes = None
        self.restype = None
        self.n_threads = None
        self.coverage_sigma = None
        self.collect_moments = None
        self.show_progress = None
        self.accumulator_cache_size = None
        self.mean_v_pointer = None
        self.cov_v_pointer = None

    def __call__(self, *args):
        self.coverage_sigma = args[8].value
        self.collect_moments = args[16].value
        self.show_progress = args[17].value
        self.accumulator_cache_size = args[18].value
        self.mean_v_pointer = args[27]
        self.cov_v_pointer = args[28]
        self.n_threads = args[-1]._obj.value
        args[19][:] = 1.0
        args[20][:] = 2.0
        args[21][:] = 3.0
        args[22][:] = 4.0
        args[23][:] = 5.0
        args[24][:] = 6.0
        args[25][:] = 7
        args[26][:] = 7.0
        args[29][:] = 0
        args[30][:] = 1.0


class _FakeGetProbabilitiesOctreeFunction:
    def __init__(self):
        self.argtypes = None
        self.restype = None
        self.n_threads = None
        self.nspatial = None
        self.max_samples_per_cell = None
        self.max_leaves_per_cell = None

    def __call__(self, *args):
        self.nspatial = args[7].value
        self.max_samples_per_cell = args[18].value
        self.max_leaves_per_cell = args[19].value
        self.n_threads = args[-1]._obj.value

        sample_stride = self.max_samples_per_cell
        leaf_stride = self.max_leaves_per_cell

        sample_spatial_index = args[24]
        velocities = args[25]
        probabilities = args[26]
        leaf_spatial_index = args[27]
        leaf_bounds = args[28]
        leaf_value_min = args[29]
        leaf_value_max = args[30]
        leaf_depth = args[31]
        leaf_sample_start = args[32]
        leaf_sample_count = args[33]
        status = args[34]
        sample_count = args[35]
        leaf_count = args[36]
        actual_sample_count = args[37]._obj
        actual_leaf_count = args[38]._obj

        sample_spatial_index[:] = -1
        velocities[:] = 0.0
        probabilities[:] = -1.0
        leaf_spatial_index[:] = -1
        leaf_bounds[:] = 0.0
        leaf_value_min[:] = 0.0
        leaf_value_max[:] = 0.0
        leaf_depth[:] = 0
        leaf_sample_start[:] = 0
        leaf_sample_count[:] = 0

        sample_spatial_index[0] = 0
        velocities[0, :] = [1.0, 2.0, 3.0]
        probabilities[0] = 0.25
        sample_spatial_index[1] = 0
        velocities[1, :] = [4.0, 5.0, 6.0]
        probabilities[1] = 0.5

        second_offset = sample_stride
        sample_spatial_index[second_offset] = 1
        velocities[second_offset, :] = [7.0, 8.0, 9.0]
        probabilities[second_offset] = -1.0

        leaf_spatial_index[0] = 0
        leaf_bounds[0, :] = [-1.0, 1.0, -2.0, 2.0, -3.0, 3.0]
        leaf_value_min[0] = 0.25
        leaf_value_max[0] = 0.5
        leaf_depth[0] = 1
        leaf_sample_start[0] = 0
        leaf_sample_count[0] = 2

        second_leaf_offset = leaf_stride
        leaf_spatial_index[second_leaf_offset] = 1
        leaf_bounds[second_leaf_offset, :] = [-4.0, 4.0, -5.0, 5.0, -6.0, 6.0]
        leaf_value_min[second_leaf_offset] = 0.0
        leaf_value_max[second_leaf_offset] = 0.0
        leaf_depth[second_leaf_offset] = 0
        leaf_sample_start[second_leaf_offset] = 0
        leaf_sample_count[second_leaf_offset] = 1

        status[:] = [0, 1]
        sample_count[:] = [2, 1]
        leaf_count[:] = [1, 1]
        actual_sample_count.value = 3
        actual_leaf_count.value = 2


class _FakeDll:
    def __init__(self):
        self.get_backtraces = _FakeBacktracesFunction()
        self.estimate_velocity_range_map = _FakeEstimateVelocityRangeMapFunction()
        self.get_probabilities_octree = _FakeGetProbabilitiesOctreeFunction()


class WrapperTypingTest(unittest.TestCase):
    def test_public_backtrace_return_annotations_match_runtime_shapes(self):
        hints = get_type_hints(wrapper.get_backtrace)

        self.assertEqual(
            hints["return"],
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        )

    def test_public_backtraces_return_annotations_match_runtime_shapes(self):
        hints = get_type_hints(wrapper.get_backtraces)

        self.assertEqual(
            hints["return"],
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        )

    def test_autorange_public_defaults_are_nonadaptive(self):
        estimate_signature = inspect.signature(wrapper.estimate_velocity_range_map)
        validate_signature = inspect.signature(wrapper.validate_and_expand_velocity_range_map)

        self.assertIs(estimate_signature.parameters["use_adaptive_dt"].default, False)
        self.assertIs(validate_signature.parameters["use_adaptive_dt"].default, False)


class GetBacktracesDllTest(unittest.TestCase):
    def test_omitted_n_threads_defaults_to_one(self):
        fake_dll = _FakeDll()

        with patch.object(wrapper.emout, "Emout", return_value=_FakeData()), \
             patch.object(
                 wrapper,
                 "create_relocated_ebvalues",
                 return_value=np.zeros((2, 2, 2, 9), dtype=np.float64),
             ), \
             patch.object(wrapper, "TempolaryInput", _FakeTemporaryInput):
            wrapper.get_backtraces_dll(
                directory="unused",
                ispec=0,
                istep=0,
                particles=[Particle([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])],
                dt=0.1,
                max_step=2,
                output_interval=1,
                use_adaptive_dt=False,
                max_probability_types=100,
                dll=fake_dll,
            )

        self.assertEqual(fake_dll.get_backtraces.n_threads, 1)


class EstimateVelocityRangeMapDllTest(unittest.TestCase):
    def test_returns_velocity_range_map_from_fortran_buffers(self):
        fake_dll = _FakeDll()

        with patch.object(wrapper.emout, "Emout", return_value=_FakeData()), \
             patch.object(
                 wrapper,
                 "create_relocated_ebvalues",
                 return_value=np.zeros((2, 2, 2, 9), dtype=np.float64),
             ), \
             patch.object(wrapper, "TempolaryInput", _FakeTemporaryInput):
            range_map = wrapper.estimate_velocity_range_map_dll(
                directory="unused",
                ispec=0,
                istep=0,
                dt=0.25,
                max_step=4,
                use_adaptive_dt=True,
                coverage_sigma=4.0,
                safety_factor=1.25,
                max_probability_types=100,
                dll=fake_dll,
                n_threads=2,
                show_progress=False,
                accumulator_cache_size=123,
            )

        self.assertEqual(fake_dll.estimate_velocity_range_map.n_threads, 2)
        self.assertEqual(fake_dll.estimate_velocity_range_map.coverage_sigma, 4.0)
        self.assertEqual(fake_dll.estimate_velocity_range_map.collect_moments, 0)
        self.assertEqual(fake_dll.estimate_velocity_range_map.show_progress, 0)
        self.assertEqual(fake_dll.estimate_velocity_range_map.accumulator_cache_size, 123)
        self.assertFalse(bool(fake_dll.estimate_velocity_range_map.mean_v_pointer))
        self.assertFalse(bool(fake_dll.estimate_velocity_range_map.cov_v_pointer))
        self.assertEqual(range_map.x_edges.tolist(), [0.0, 1.0])
        np.testing.assert_allclose(range_map.vx_min, [[[1.0]]])
        np.testing.assert_allclose(range_map.vx_max, [[[2.0]]])
        np.testing.assert_array_equal(range_map.count, [[[7]]])
        self.assertEqual(range_map.mean_v.shape, (1, 1, 1, 3))
        self.assertTrue(np.isnan(range_map.mean_v).all())
        np.testing.assert_allclose(range_map.confidence, [[[1.0]]])
        self.assertEqual(range_map.directory, _FakeData.directory)
        self.assertEqual(range_map.metadata["directory"], str(_FakeData.directory))
        self.assertEqual(range_map.metadata["ispec"], 0)
        self.assertEqual(range_map.metadata["istep"], 0)


class GetProbabilitiesOctreeDllTest(unittest.TestCase):
    def test_phase_grid_call_returns_compacted_octree_result(self):
        fake_dll = _FakeDll()
        phase_grid = PhaseGrid(
            x=(0.0, 1.0, 2),
            y=0.0,
            z=0.0,
            vx=(-1.0, 1.0, 2),
            vy=(-2.0, 2.0, 2),
            vz=(-3.0, 3.0, 2),
        )

        with patch.object(wrapper.emout, "Emout", return_value=_FakeData()), \
             patch.object(
                 wrapper,
                 "create_relocated_ebvalues",
                 return_value=np.zeros((2, 2, 2, 9), dtype=np.float64),
             ), \
             patch.object(wrapper, "TempolaryInput", _FakeTemporaryInput):
            result = wrapper.get_probabilities_octree_dll(
                directory="unused",
                ispec=0,
                istep=0,
                phase_grid=phase_grid,
                position=None,
                velocity_bounds=None,
                dt=0.25,
                max_step=4,
                use_adaptive_dt=False,
                max_probability_types=100,
                scout_bins=(3, 3, 3),
                max_depth=2,
                max_samples_per_cell=3,
                max_leaves_per_cell=2,
                refine_threshold_rel=1e-4,
                edge_threshold_rel=1e-4,
                expand_factor=1.5,
                max_expansions=1,
                dll=fake_dll,
                n_threads=2,
            )

        self.assertEqual(fake_dll.get_probabilities_octree.n_threads, 2)
        self.assertEqual(fake_dll.get_probabilities_octree.nspatial, 2)
        self.assertEqual(fake_dll.get_probabilities_octree.max_samples_per_cell, 3)
        np.testing.assert_allclose(result.spatial_points, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        np.testing.assert_array_equal(result.spatial_index, [0, 0, 1])
        np.testing.assert_allclose(result.velocities, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        self.assertTrue(np.isnan(result.probabilities[-1]))
        np.testing.assert_array_equal(result.leaf_spatial_index, [0, 1])
        np.testing.assert_array_equal(result.leaf_sample_start, [0, 2])
        self.assertEqual(result.leaf_sample_slice(1), slice(2, 3))
        np.testing.assert_array_equal(result.sample_count, [2, 1])
        np.testing.assert_array_equal(result.leaf_count, [1, 1])
        self.assertEqual(result.metadata["actual_sample_count"], 3)
        self.assertEqual(result.metadata["actual_leaf_count"], 2)

    def test_rejects_ambiguous_octree_input_forms(self):
        with self.assertRaises(ValueError):
            wrapper._prepare_octree_inputs(
                phase_grid=PhaseGrid(0.0, 0.0, 0.0, (-1.0, 1.0, 2), (-1.0, 1.0, 2), (-1.0, 1.0, 2)),
                position=[0.0, 0.0, 0.0],
                velocity_bounds=((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)),
            )

    def test_accepts_per_spatial_octree_velocity_bounds(self):
        spatial_points, velocity_bounds = wrapper._prepare_octree_inputs(
            phase_grid=None,
            position=[[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]],
            velocity_bounds=[
                [-1.0, 1.0, -2.0, 2.0, -3.0, 3.0],
                [-4.0, 4.0, -5.0, 5.0, -6.0, 6.0],
            ],
        )

        np.testing.assert_allclose(spatial_points, [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
        np.testing.assert_allclose(
            velocity_bounds,
            [
                [-1.0, 1.0, -2.0, 2.0, -3.0, 3.0],
                [-4.0, 4.0, -5.0, 5.0, -6.0, 6.0],
            ],
        )

    def test_broadcasts_single_octree_velocity_bounds_to_all_positions(self):
        spatial_points, velocity_bounds = wrapper._prepare_octree_inputs(
            phase_grid=None,
            position=[[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]],
            velocity_bounds=[[-1.0, 1.0], [-2.0, 2.0], [-3.0, 3.0]],
        )

        np.testing.assert_allclose(spatial_points, [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
        np.testing.assert_allclose(
            velocity_bounds,
            [
                [-1.0, 1.0, -2.0, 2.0, -3.0, 3.0],
                [-1.0, 1.0, -2.0, 2.0, -3.0, 3.0],
            ],
        )

    def test_accepts_per_spatial_octree_bounds_as_axis_pairs(self):
        _, velocity_bounds = wrapper._prepare_octree_inputs(
            phase_grid=None,
            position=[[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]],
            velocity_bounds=[
                [[-1.0, 1.0], [-2.0, 2.0], [-3.0, 3.0]],
                [[-4.0, 4.0], [-5.0, 5.0], [-6.0, 6.0]],
            ],
        )

        np.testing.assert_allclose(
            velocity_bounds,
            [
                [-1.0, 1.0, -2.0, 2.0, -3.0, 3.0],
                [-4.0, 4.0, -5.0, 5.0, -6.0, 6.0],
            ],
        )


class VelocityRangeValidationTest(unittest.TestCase):
    def test_edge_probability_mask_flags_only_cells_with_edge_signal(self):
        prob_grid = np.zeros((1, 1, 2, 3, 3, 3), dtype=np.float64)
        prob_grid[0, 0, 0, 1, 1, 1] = 1.0
        prob_grid[0, 0, 0, 0, 1, 1] = 0.01
        prob_grid[0, 0, 1, 1, 1, 1] = 1.0
        prob_grid[0, 0, 1, 2, 1, 1] = 0.2

        mask = wrapper._edge_probability_mask(prob_grid, edge_threshold=0.1)

        np.testing.assert_array_equal(mask, [[[False, True]]])


if __name__ == "__main__":
    unittest.main()
