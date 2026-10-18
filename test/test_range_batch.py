"""Tier-2 bundles, Python-free executable parity, persistence and resume."""
import copy
import json
import os
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

import f90nml
import numpy as np

from test_batch import batch_case, EXECUTABLE
from vdsolverf.core import VelocityRangeMap
from vdsolverf.emses import wrapper
from vdsolverf.io import write_range_batch, read_range_batch_results

# Use the installed shared library when the source tree has none (sdist CI).
from test_native_regressions import LIBRARY_PATH
from vdsolverf.emses import wrapper as _wrapper
_wrapper._DEFAULT_LIBRARY_PATHS["linux"] = LIBRARY_PATH


FIELDS = ("vx_min", "vx_max", "vy_min", "vy_max", "vz_min", "vz_max",
          "weight_sum", "confidence", "count", "status")
SETTINGS = dict(dt=0.7, max_step=12, coverage_sigma=2.3, safety_factor=1.2,
                source_samples_per_cell=1, minimum_count=3, accumulator_cache_size=7, n_threads=1)


class RangeBundleTest(unittest.TestCase):
    def compare(self, actual, expected):
        for name in FIELDS:
            np.testing.assert_array_equal(getattr(actual, name), getattr(expected, name), err_msg=name)

    def test_bundle_and_live_estimation_and_validation_are_identical(self):
        with batch_case() as (_, directory, bundle, _, _, _):
            # Use the same noncubic nonzero relocated fields for the live path,
            # rather than the base fixture's smaller zero-field stub.
            with patch.object(wrapper, "create_relocated_ebvalues", return_value=bundle[0]):
                live = wrapper.estimate_velocity_range_map(directory, 0, 0, show_progress=False,
                                                           collect_moments=True, **SETTINGS)
                validated = wrapper.validate_and_expand_velocity_range_map(
                    copy.deepcopy(live), directory, 0, 0, dt=0.7, max_step=12,
                    coarse_bins=(3, 3, 3), edge_threshold=0.1, max_iter=1, n_threads=1)
            with patch.object(wrapper.emout, "Emout", side_effect=AssertionError("must use bundle")), \
                 patch.object(wrapper, "create_relocated_ebvalues", side_effect=AssertionError("must use bundle")), \
                 patch.object(wrapper, "TempolaryInput", side_effect=AssertionError("must use bundle")):
                for source in (bundle, bundle.namelist_path.parent):
                    mapped = wrapper.estimate_velocity_range_map(directory, 0, 0, field_bundle=source,
                                                                 collect_moments=True, show_progress=False, **SETTINGS)
                    self.compare(mapped, live)
                    np.testing.assert_array_equal(mapped.mean_v, live.mean_v)
                    np.testing.assert_array_equal(mapped.cov_v, live.cov_v)
                    checked = wrapper.validate_and_expand_velocity_range_map(
                        mapped, directory, 0, 0, dt=0.7, max_step=12, coarse_bins=(3, 3, 3),
                        edge_threshold=0.1, max_iter=1, n_threads=1, field_bundle=source)
                    self.compare(checked, validated)

    def test_consistency_checks_and_writer_settings(self):
        with batch_case() as (_, directory, bundle, _, _, _):
            live = wrapper.estimate_velocity_range_map(directory, 0, 0, field_bundle=bundle,
                                                       show_progress=False, **SETTINGS)
            for options in (dict(directory=directory / "other", ispec=0, istep=0),
                            dict(directory=directory, ispec=0, istep=1),
                            dict(directory=directory, ispec=1, istep=0),
                            dict(directory=directory, ispec=0, istep=0, use_electric_field=False),
                            dict(directory=directory, ispec=0, istep=0, use_magnetic_field=False)):
                with self.subTest(options=options):
                    with self.assertRaises(ValueError):
                        wrapper.estimate_velocity_range_map(field_bundle=bundle, **options, **SETTINGS)
                    with self.assertRaises(ValueError):
                        wrapper.validate_and_expand_velocity_range_map(live, field_bundle=bundle,
                                                                      dt=0.7, max_step=12, max_iter=0, **options)
            for options in (dict(coverage_sigma=np.nan), dict(safety_factor=0.9), dict(n_threads=0),
                            dict(minimum_count=0), dict(source_samples_per_cell=0), dict(dt=0),
                            dict(accumulator_cache_size=0), dict(max_step=1.5), dict(typo=1),
                            dict(coverage_mode="bad"), dict(coverage_mode="relative_density", eps_rel=np.nan),
                            dict(use_electric_field=False)):
                with self.subTest(options=options), self.assertRaises(ValueError):
                    write_range_batch(directory / "bad", bundle, 0, **options)
            path = write_range_batch(directory / "ranges", bundle, 0,
                                     coverage_mode="relative_density", eps_rel=1e-6, n_threads=None)
            config = f90nml.read(path)["vdsolverf_batch"]
            self.assertEqual(config["mode"], "velocity_range")
            self.assertEqual(config["dt"], 0.25)
            self.assertAlmostEqual(config["coverage_sigma"], np.sqrt(-2*np.log(1e-6)))
            self.assertNotIn("coverage_mode", config)
            self.assertNotIn("n_threads", config)
            self.assertEqual({p.name for p in path.parent.iterdir()}, {"batch.nml"})
            with self.assertRaisesRegex(ValueError, "incomplete"):
                read_range_batch_results(path)


@unittest.skipUnless(EXECUTABLE.is_file(), "requires make batch")
class RangeBatchExecutableTest(unittest.TestCase):
    def run_batch(self, path, *args, env=None):
        result = subprocess.run([str(EXECUTABLE.resolve()), *args, str(path)], cwd=path.parent,
                                text=True, capture_output=True, timeout=90, env=env)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout

    def test_serial_matches_api_bit_for_bit_and_save_roundtrip(self):
        with batch_case() as (_, directory, bundle, _, _, _):
            reference = wrapper.estimate_velocity_range_map(directory, 0, 0, field_bundle=bundle,
                                                            show_progress=False, **SETTINGS)
            path = write_range_batch(directory / "ranges", bundle, 0, **SETTINGS)
            self.assertIn("computed", self.run_batch(path))
            result = read_range_batch_results(path)
            for name in FIELDS:
                actual, expected = getattr(result, name), getattr(reference, name)
                np.testing.assert_array_equal(actual, expected, err_msg=name)
                # Compare raw representations, including NaNs, with GNU builds.
                np.testing.assert_array_equal(actual.view(np.uint8), expected.view(np.uint8), err_msg=name)
                self.assertTrue(actual.flags.c_contiguous)
                self.assertEqual(actual.dtype, np.int32 if name in ("count", "status") else np.float64)
            self.assertEqual(result.cell_shape, (4, 2, 3))
            self.assertEqual(result.metadata["istep"], 0)
            self.assertEqual(result.metadata["ispec"], 0)
            self.assertEqual(result.metadata["n_threads"], 1)
            self.assertTrue(result.metadata["use_electric_field"])
            self.assertTrue(result.metadata["use_magnetic_field"])
            self.assertEqual(sum(result.metadata["summary"]["status_counts"].values()), 24)
            loaded = VelocityRangeMap.load(result.save())
            for name in FIELDS:
                np.testing.assert_array_equal(getattr(result, name), getattr(loaded, name))
            self.assertFalse(loaded._mean_v_stored)
            self.assertFalse(loaded._cov_v_stored)
            bad = path.parent / "out/velocity_range/count.npy"
            np.save(bad, np.zeros((4, 2, 3), dtype=np.int64))
            with self.assertRaisesRegex(ValueError, "C-order <i4"):
                read_range_batch_results(path)

    def test_resume_checks_all_input_contents_and_missing_outputs(self):
        with batch_case() as (_, directory, bundle, _, _, _):
            path = write_range_batch(directory / "ranges", bundle, 0, **SETTINGS)
            self.run_batch(path)
            output = path.parent / "out/velocity_range"
            before = {p.name: p.stat().st_mtime_ns for p in output.iterdir()}
            self.assertIn("skipped", self.run_batch(path, "--resume"))
            self.assertEqual(before, {p.name: p.stat().st_mtime_ns for p in output.iterdir()})
            previous = json.loads((output / "summary.json").read_text())["input_fingerprint"]
            for changed in (path, bundle.namelist_path, bundle.namelist_path.parent / "bundle.json",
                            bundle.namelist_path.parent / "ebvalues_ispec0.npy"):
                stamp = changed.stat()
                if changed.suffix == ".npy":
                    values = np.load(changed)
                    values.flat[-1] += 1e-5
                    np.save(changed, values)
                else:
                    changed.write_text(changed.read_text() + "\n")
                os.utime(changed, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
                self.assertIn("computed", self.run_batch(path))
                current = json.loads((output / "summary.json").read_text())["input_fingerprint"]
                self.assertNotEqual(previous, current)
                previous = current
            (output / "count.npy").unlink()
            self.assertIn("computed", self.run_batch(path))
            self.assertIn("computed", self.run_batch(path, "--no-resume"))
            self.assertIn("skipped", self.run_batch(path))

    def test_runtime_threads_and_bundle_field_flags_are_inherited(self):
        with batch_case() as (_, directory, bundle, _, _, _):
            metadata_path = bundle.namelist_path.parent / "bundle.json"
            metadata = json.loads(metadata_path.read_text())
            metadata.update(use_electric_field=False, use_magnetic_field=False)
            metadata_path.write_text(json.dumps(metadata))
            path = write_range_batch(directory / "ranges", metadata_path.parent, 0,
                                     **dict(SETTINGS, n_threads=None))
            for threads in (3, 4):
                self.assertIn("computed", self.run_batch(path, env=dict(os.environ, OMP_NUM_THREADS=str(threads))))
                result = read_range_batch_results(path)
                self.assertEqual(result.metadata["n_threads"], threads)
                self.assertFalse(result.metadata["use_electric_field"])
                self.assertFalse(result.metadata["use_magnetic_field"])


if __name__ == "__main__":
    unittest.main()
