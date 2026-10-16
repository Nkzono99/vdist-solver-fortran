"""Batch helpers, real serial executable, checkpoint and input-format regressions."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

import f90nml
import numpy as np

from test_dust_currents import maxwellian_input
from vdsolverf.core import DustCurrentResult
from vdsolverf.emses import wrapper
from vdsolverf.io import read_batch_results, write_batch, write_field_bundle, write_batch_results


ROOT = Path(__file__).resolve().parents[1]
EXECUTABLE = Path(os.environ.get("VDSOLVERF_BATCH_EXECUTABLE", ROOT / "build/batch-serial/bin/vdsolverf-batch"))
FIELDS = ("flux", "flux_error", "status", "sample_count", "terminal_box_count")


@contextmanager
def batch_case(**options):
    with maxwellian_input() as (directory, data, _):
        # Noncubic and nonzero, with distinct components at each grid point:
        # the C -> Fortran axis reversal must match ctypes exactly.
        data.inp.nml["tmgrid"].update(nx=3, ny=2, nz=4)
        fields = np.arange(5 * 3 * 4 * 9, dtype=np.float64).reshape(5, 3, 4, 9) * 1e-5
        with patch.object(wrapper, "create_relocated_ebvalues", return_value=fields):
            bundle = write_field_bundle(directory, 0, [0], directory / "bundle")
        targets = np.array([[0.7 + 0.1*i, 0.8, 1.7, -4., 4., -4., 4., -4., 4.] for i in range(7)])
        targets[-1, 4] = -5.  # preserve native invalid-box status 5
        psi = np.array([-1., 0., 0.25, 1.])
        solver = dict(dt=0.7, max_step=32, scout_bins=[3, 4, 3], max_depth=1,
                      max_samples_per_cell=4096, max_leaves_per_cell=64,
                      max_expansions=0, chunk_size=2, n_threads=2)
        solver.update(options)
        path = write_batch(directory / "batch", bundle, 0, targets, psi, **solver)
        yield path, directory, bundle, targets, psi, solver


def write_mock_results(path, targets, psi):
    arrays = dict(flux=np.arange(len(targets)*len(psi), dtype="<f8").reshape(len(targets), len(psi)),
                  flux_error=np.ones((len(targets), len(psi)), dtype="<f8"),
                  status=np.arange(len(targets), dtype="<i4") % 6,
                  sample_count=np.arange(len(targets), dtype="<i8") + 2**40,
                  terminal_box_count=np.arange(len(targets), dtype="<i8") + 2**35)
    write_batch_results(path, DustCurrentResult(targets[:, :3], psi, **arrays))
    return arrays


class BatchHelperTest(unittest.TestCase):
    def test_write_and_collect_roundtrip_and_int64_diagnostics(self):
        with batch_case() as (path, _, bundle, targets, psi, _):
            config = f90nml.read(path)["vdsolverf_batch"]
            self.assertEqual(config["bundle_dir"], "../bundle")
            self.assertEqual(config["ispec"], 0)
            self.assertEqual(config["chunk_size"], 2)
            np.testing.assert_array_equal(np.load(path.parent / "targets.npy"), targets)
            arrays = write_mock_results(path, targets, psi)
            result = read_batch_results(path.parent)
            self.assertIsInstance(result, DustCurrentResult)
            np.testing.assert_array_equal(result.positions, targets[:, :3])
            np.testing.assert_array_equal(result.psi_grid, psi)
            for name in FIELDS:
                np.testing.assert_array_equal(getattr(result, name), arrays[name])
            self.assertEqual(result.sample_count.dtype, np.int64)
            self.assertEqual(result.terminal_box_count.dtype, np.int64)
            self.assertEqual(result.status.dtype, np.int32)
            # dt defaults from the bundle and n_threads stays runtime-dependent.
            other = write_batch(path.parent / "other", bundle, 0, targets, psi, n_threads=None)
            settings = f90nml.read(other)["vdsolverf_batch"]
            self.assertEqual(settings["dt"], 0.1)
            self.assertNotIn("n_threads", settings)

    def test_collection_rejects_missing_markers_wrong_dtypes_and_shapes(self):
        with batch_case() as (path, _, _, targets, psi, _):
            write_mock_results(path, targets, psi)
            marker = path.parent / "out/chunk_000000002.done"
            marker.unlink()
            with self.assertRaisesRegex(ValueError, "incomplete"):
                read_batch_results(path)
            marker.touch()
            bad = path.parent / "out/chunk_000000002_status.npy"
            np.save(bad, np.zeros(2, dtype=np.int64))
            with self.assertRaisesRegex(ValueError, "C-order <i4"):
                read_batch_results(path)
            np.save(bad, np.zeros(1, dtype=np.int32))
            with self.assertRaisesRegex(ValueError, "invalid chunk"):
                read_batch_results(path)

    def test_prepare_rejects_invalid_settings_and_field_mismatches(self):
        with batch_case() as (path, _, bundle, targets, psi, _):
            for settings in ({"chunk_size": 0}, {"max_depth": -1}, {"scout_bins": [1, 2, 3]},
                             {"dt": np.nan}, {"refine_threshold_rel": np.inf}, {"n_threads": 0},
                             {"typo": 1}, {"use_electric_field": False}, {"use_magnetic_field": False}):
                with self.subTest(settings=settings), self.assertRaises(ValueError):
                    write_batch(path.parent / "invalid", bundle, 0, targets, psi, **settings)
            with self.assertRaises(ValueError):
                write_batch(path.parent / "invalid", bundle, 1, targets, psi)
            for bad in ([], [[0]*8], np.full((1, 9), np.nan)):
                with self.assertRaises(ValueError):
                    write_batch(path.parent / "invalid", bundle, 0, bad, psi)
            for bad in ([], [np.nan], [[0.]]):
                with self.assertRaises(ValueError):
                    write_batch(path.parent / "invalid", bundle, 0, targets, bad)


@unittest.skipUnless(EXECUTABLE.is_file(), "serial vdsolverf-batch has not been built (make batch)")
class BatchExecutableTest(unittest.TestCase):
    def run_batch(self, path, *args, success=True):
        result = subprocess.run([str(EXECUTABLE.resolve()), *args, str(path)], cwd=path.parent.parent,
                                text=True, capture_output=True, timeout=90)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_ctypes_equality_and_resume_recovery(self):
        with batch_case() as (path, directory, bundle, targets, psi, solver):
            library_solver = {key: value for key, value in solver.items() if key != "chunk_size"}
            reference = wrapper.get_dust_currents(directory, 0, 0, targets[:, :3], targets[:, 3:], psi,
                                                  field_bundle=bundle, **library_solver)
            self.run_batch(path)
            result = read_batch_results(path)
            for name in FIELDS:
                if name in ("flux", "flux_error"):
                    np.testing.assert_allclose(getattr(result, name), getattr(reference, name), rtol=2e-14, atol=1e-14)
                else:
                    np.testing.assert_array_equal(getattr(result, name), getattr(reference, name))
            self.assertEqual(result.status[-1], 5)
            self.assertEqual(result.metadata["summary"]["chunks_computed"], 4)
            self.assertEqual(result.metadata["summary"]["status_counts"]["5"], 1)
            output = path.parent / "out"
            before = {p.name: p.stat().st_mtime_ns for p in output.glob("chunk_*")}
            self.run_batch(path, "--resume")
            resumed = read_batch_results(path)
            self.assertEqual(resumed.metadata["summary"]["chunks_computed"], 0)
            self.assertEqual(resumed.metadata["summary"]["chunks_skipped"], 4)
            self.assertEqual(before, {p.name: p.stat().st_mtime_ns for p in output.glob("chunk_*")})
            (output / "chunk_000000002.done").unlink()
            np.save(output / "chunk_000000002_flux.npy", np.array([-99.]))
            self.run_batch(path)
            recovered = read_batch_results(path)
            self.assertEqual(recovered.metadata["summary"]["chunks_computed"], 1)
            self.assertEqual(recovered.metadata["summary"]["chunks_skipped"], 3)
            np.testing.assert_array_equal(recovered.flux, result.flux)
            self.run_batch(path, "--no-resume")
            self.assertEqual(read_batch_results(path).metadata["summary"]["chunks_computed"], 4)

    def test_npy_v2_inputs_and_fields(self):
        with batch_case() as (path, _, bundle, targets, psi, _):
            self.run_batch(path)
            reference = read_batch_results(path)
            for name in (path.parent / "targets.npy", path.parent / "psi.npy",
                         bundle.namelist_path.parent / "ebvalues_ispec0.npy"):
                values = np.load(name)
                with name.open("wb") as stream:
                    np.lib.format.write_array(stream, values, version=(2, 0), allow_pickle=False)
            self.run_batch(path, "--no-resume")
            result = read_batch_results(path)
            for name in FIELDS:
                np.testing.assert_array_equal(getattr(result, name), getattr(reference, name))

    def test_default_threads_follow_execution_environment(self):
        with batch_case() as (path, _, bundle, targets, psi, solver):
            solver["n_threads"] = None
            write_batch(path.parent, bundle, 0, targets, psi, **solver)
            result = subprocess.run([str(EXECUTABLE.resolve()), str(path)], text=True, capture_output=True,
                                    timeout=90, env=dict(os.environ, OMP_NUM_THREADS="3"))
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(read_batch_results(path).metadata["summary"]["threads"], 3)

    def test_rejects_unsupported_npy_and_bundle_flags(self):
        with batch_case() as (path, _, bundle, targets, _, _):
            target_path = path.parent / "targets.npy"
            for values, message in ((targets.astype(np.float32), "dtype"),
                                    (targets.astype(">f8"), "little-endian"),
                                    (np.asfortranarray(targets), "Fortran-order"),
                                    (targets[:, :8], "shape")):
                with self.subTest(message=message):
                    np.save(target_path, values)
                    result = self.run_batch(path, success=False)
                    self.assertIn(message, result.stderr)
            np.save(target_path, targets)
            raw = target_path.read_bytes()
            target_path.write_bytes(raw[:-1])
            self.assertIn("payload length", self.run_batch(path, success=False).stderr)
            target_path.write_bytes(raw)
            metadata_path = bundle.namelist_path.parent / "bundle.json"
            metadata = json.loads(metadata_path.read_text())
            metadata["use_electric_field"] = False
            metadata_path.write_text(json.dumps(metadata))
            self.assertIn("use_electric_field", self.run_batch(path, success=False).stderr)


if __name__ == "__main__":
    unittest.main()
