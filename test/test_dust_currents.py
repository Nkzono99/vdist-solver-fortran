"""Real ctypes and bundle regressions on an isotropic open-boundary Maxwellian."""
from contextlib import ExitStack, contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import f90nml
import numpy as np

from vdsolverf.core import DustCurrentResult
from vdsolverf.emses import get_dust_currents, get_probabilities_octree, wrapper
from vdsolverf.io import FieldBundle, read_field_bundle, write_field_bundle

# Use the installed shared library when the source tree has none (sdist CI).
from test_native_regressions import LIBRARY_PATH
from vdsolverf.emses import wrapper as _wrapper
_wrapper._DEFAULT_LIBRARY_PATHS["linux"] = LIBRARY_PATH


@contextmanager
def maxwellian_input():
    with TemporaryDirectory() as directory, ExitStack() as stack:
        directory = Path(directory)
        nml = f90nml.Namelist({
            'esorem': {'emflag': 0},
            'plasma': {'wp': [0.] * 10, 'wc': 0., 'phiz': 0., 'phixy': 0.},
            'tmgrid': {'nx': 2, 'ny': 2, 'nz': 2, 'dt': 0.1},
            'system': {'nspec': 1, 'npbnd': [[2, 2, 2]] * 10},
            'intp': {'qm': [1.] * 10, 'peth': [1.], 'path': [1.],
                     'spa': [0.] * 10, 'spe': [0.] * 10, 'speth': [0.] * 10,
                     'vdri': [0.] * 10, 'vdthz': [0.] * 10, 'vdthxy': [0.] * 10},
            'ptcond': {'boundary_type': 'none'},
            'emissn': {'nepl': [0] * 10, 'nflag_emit': [0] * 10},
        })
        class Input:
            def __init__(self):
                self.nml = nml
            def __contains__(self, key):
                return any(key in group for group in nml.values())
            def __getattr__(self, key):
                for group in nml.values():
                    if key in group:
                        return group[key]
                raise AttributeError(key)
        data = SimpleNamespace(directory=directory, inp=Input())
        fields = np.zeros((3, 3, 3, 9), dtype=np.float64)
        stack.enter_context(patch.object(wrapper.emout, 'Emout', return_value=data))
        stack.enter_context(patch.object(wrapper, 'create_relocated_ebvalues', return_value=fields))
        yield directory, data, fields


class DustCurrentTest(unittest.TestCase):
    def dust(self, directory, **kwargs):
        options = dict(positions=[[1., 1., 1.], [0.75, 1., 1.]],
                       velocity_bounds=[-4.5, 4.5] * 3, psi_grid=[-1., -0.25, 0., 0.25, 1.],
                       dt=2., max_step=32, n_threads=1)
        options.update(kwargs)
        return get_dust_currents(directory, 0, 0, **options)

    def octree(self, directory, **kwargs):
        options = dict(position=[[1., 1., 1.], [0.75, 1., 1.]],
                       velocity_bounds=[-4.5, 4.5] * 3, dt=2., max_step=32,
                       max_depth=2, n_threads=1)
        options.update(kwargs)
        return get_probabilities_octree(directory, 0, 0, **options)

    def test_native_maxwellian_default_accuracy_and_parallel_equality(self):
        with maxwellian_input() as (directory, _, _):
            result = self.dust(directory)
            self.assertIsInstance(result, DustCurrentResult)
            psi = result.psi_grid
            expected = np.sqrt(8 / np.pi) * np.where(psi > 0, np.exp(-psi), 1 - psi)
            np.testing.assert_allclose(result.flux, np.broadcast_to(expected, result.flux.shape), rtol=1e-2)
            self.assertTrue(np.all(np.isfinite(result.flux_error)))
            self.assertTrue(np.all(result.flux_error >= 0))
            parallel = self.dust(directory, n_threads=2)
            for name in ('flux', 'flux_error', 'status', 'sample_count', 'terminal_box_count'):
                np.testing.assert_array_equal(getattr(result, name), getattr(parallel, name))

    def test_bundle_roundtrip_and_exact_native_equality(self):
        with maxwellian_input() as (directory, data, fields):
            from vdsolverf.emses.tmpolary_input import TempolaryInput
            with TempolaryInput(data) as inp:
                expected_namelist = inp.tmppath.read_bytes()
            bundle = write_field_bundle(directory, 0, [0], directory / 'bundle')
            self.assertIsInstance(bundle, FieldBundle)
            self.assertEqual(bundle.namelist_path.read_bytes(), expected_namelist)
            restored = read_field_bundle(directory / 'bundle' / 'bundle.json')
            self.assertEqual(restored._arrays, {})
            self.assertIsInstance(restored[0], np.memmap)
            np.testing.assert_array_equal(restored[0], fields)
            self.assertFalse(restored[0].flags.writeable)
            with (directory / 'bundle' / 'ebvalues_ispec0.npy').open('rb') as stream:
                self.assertEqual(np.lib.format.read_magic(stream), (1, 0))
                shape, fortran_order, dtype = np.lib.format.read_array_header_1_0(stream)
            self.assertEqual(shape, fields.shape)
            self.assertFalse(fortran_order)
            self.assertEqual(dtype.str, '<f8')
            direct = self.dust(directory)
            direct_octree = self.octree(directory)
            with patch.object(wrapper.emout, 'Emout', side_effect=AssertionError('must not read EMSES')), \
                 patch.object(wrapper, 'create_relocated_ebvalues', side_effect=AssertionError('must not relocate')):
                bundled = self.dust(directory, field_bundle=restored)
                bundled_octree = self.octree(directory, field_bundle=directory / 'bundle')
            for name in ('flux', 'flux_error', 'status', 'sample_count', 'terminal_box_count'):
                np.testing.assert_array_equal(getattr(direct, name), getattr(bundled, name))
            for name in ('velocities', 'probabilities', 'leaf_bounds', 'leaf_is_terminal', 'status', 'sample_count'):
                np.testing.assert_array_equal(getattr(direct_octree, name), getattr(bundled_octree, name))

    def test_terminal_boxes_partition_expanded_root(self):
        with maxwellian_input() as (directory, _, _):
            result = self.octree(directory, velocity_bounds=[-3., 3.] * 3)
            self.assertEqual(result.leaf_is_terminal.dtype, np.dtype(bool))
            for spatial in range(2):
                mask = result.leaf_mask(spatial)
                root = result.leaf_bounds[mask & (result.leaf_depth == 0)][0]
                boxes = result.leaf_bounds[mask & result.leaf_is_terminal]
                self.assertAlmostEqual(np.prod(boxes[:, 1::2] - boxes[:, ::2], axis=1).sum(),
                                       np.prod(root[1::2] - root[::2]))
                self.assertEqual(result.status[spatial], 3)
                np.testing.assert_array_equal(root, [-4.5, 4.5] * 3)

    def test_capacity_keeps_complete_current_domain(self):
        with maxwellian_input() as (directory, _, _):
            root = self.dust(directory, scout_bins=(9, 9, 9), max_depth=0)
            # Neither capacity can hold eight children: retain the sampled root.
            for budget, code in ((dict(max_samples_per_cell=730), 2), (dict(max_leaves_per_cell=8), 4)):
                limited = self.dust(directory, **budget)
                np.testing.assert_array_equal(limited.flux, root.flux)
                np.testing.assert_array_equal(limited.terminal_box_count, [1, 1])
                np.testing.assert_array_equal(limited.status, [code, code])
            empty = self.dust(directory, max_samples_per_cell=8)
            np.testing.assert_array_equal(empty.status, [2, 2])
            np.testing.assert_array_equal(empty.flux, np.zeros((2, 5)))
            np.testing.assert_array_equal(empty.sample_count, [0, 0])

    def test_matches_terminal_box_trapezoids_and_sample_range_error(self):
        with maxwellian_input() as (directory, _, _):
            # Complete depth-two trees use the same boxes regardless of schedule.
            samples = self.octree(directory)
            moments = self.dust(directory, max_depth=2)
            expected = np.zeros_like(moments.flux)
            expected_error = np.zeros_like(moments.flux_error)
            for box in np.flatnonzero(samples.leaf_is_terminal):
                point = samples.leaf_spatial_index[box]
                selected = samples.leaf_sample_slice(box)
                speed = np.linalg.norm(samples.velocities[selected], axis=1)
                pdf = np.nan_to_num(samples.probabilities[selected])
                kernel = np.zeros((len(speed), len(moments.psi_grid)))
                moving = speed > 0
                kernel[moving] = np.maximum(0., speed[moving, None] -
                    2 * moments.psi_grid[None, :] / speed[moving, None])
                integrand = pdf[:, None] * kernel
                bounds = samples.leaf_bounds[box]
                volume = np.prod(bounds[1::2] - bounds[::2])
                w = np.array([0.5, 1., 0.5])
                weights = (w[:, None, None] * w[None, :, None] * w[None, None, :]).ravel() / 8
                expected[point] += volume * np.sum(integrand * weights[:, None], axis=0)
                expected_error[point] += volume * np.ptp(integrand, axis=0)
            np.testing.assert_allclose(moments.flux, expected, rtol=2e-14, atol=1e-15)
            np.testing.assert_allclose(moments.flux_error, expected_error, rtol=2e-14, atol=1e-15)
            np.testing.assert_array_equal(moments.terminal_box_count, [64, 64])

    def test_invalid_bounds_and_no_boundary_signal(self):
        with maxwellian_input() as (directory, _, _):
            for bounds in ([1., -1.] * 3, [-np.inf, np.inf] * 3):
                invalid = self.dust(directory, velocity_bounds=bounds)
                np.testing.assert_array_equal(invalid.status, [5, 5])
                np.testing.assert_array_equal(invalid.flux, np.zeros((2, 5)))
            absent = self.dust(directory, velocity_bounds=[-0.001, 0.001] * 3,
                               max_step=1, max_expansions=0)
            np.testing.assert_array_equal(absent.status, [1, 1])
            np.testing.assert_array_equal(absent.flux, np.zeros((2, 5)))

    def test_terminal_flags_preserve_original_octree_samples(self):
        with maxwellian_input() as (directory, _, _):
            native = self.octree(directory)
            dll = wrapper._load_dll('auto', None)
            # Hide only the new symbol to invoke the original, unchanged C ABI.
            legacy_dll = SimpleNamespace(get_probabilities_octree=dll.get_probabilities_octree)
            with patch.object(wrapper, '_load_dll', return_value=legacy_dll):
                legacy = self.octree(directory)
            for name in ('velocities', 'probabilities', 'leaf_bounds', 'leaf_is_terminal',
                         'leaf_depth', 'leaf_sample_start', 'leaf_sample_count', 'sample_count', 'status'):
                np.testing.assert_array_equal(getattr(native, name), getattr(legacy, name))
            # The old API may split a parent before its children are sampled.
            truncated = self.octree(directory, max_samples_per_cell=730)
            self.assertFalse(truncated.leaf_is_terminal[0])
            np.testing.assert_array_equal(truncated.status, [2, 2])

    def test_bundle_rejects_mismatched_request_and_bad_array(self):
        with maxwellian_input() as (directory, _, _):
            bundle = write_field_bundle(directory, 0, [0], directory / 'bundle')
            with self.assertRaises(ValueError):
                self.dust(directory, field_bundle=bundle, use_electric_field=False)
            with self.assertRaises(ValueError):
                wrapper.get_dust_currents(directory, 0, 1, [1., 1., 1.], [-4.5, 4.5]*3, [0.], 2., 32,
                                          field_bundle=bundle)
            with self.assertRaises(ValueError):
                bundle[1]
            np.save(directory / 'bundle' / 'ebvalues_ispec0.npy', np.zeros((1,), dtype=np.float32))
            with self.assertRaises(ValueError):
                read_field_bundle(directory / 'bundle')[0]

    def test_rejects_nonfinite_psi_before_native_call(self):
        with self.assertRaises(ValueError):
            self.dust('.', psi_grid=[np.nan])


if __name__ == '__main__':
    unittest.main()
