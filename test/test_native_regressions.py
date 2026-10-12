"""Exercise the real shared library with small synthetic EMSES inputs."""

from contextlib import ExitStack, contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from vdsolverf.core import Particle
from vdsolverf.emses import wrapper


LIBRARY_PATH = wrapper.VDIST_SOLVER_FORTRAN_LIBRARY_PATH_LINUX
if not LIBRARY_PATH.is_file():
    # Source-distribution CI installs the library in site-packages.
    try:
        from importlib.metadata import distribution

        LIBRARY_PATH = Path(
            distribution("vdist-solver-fortran").locate_file(
                "vdsolverf/libvdist-solver-fortran.so"
            )
        )
    except (ImportError, LookupError):
        pass


@contextmanager
def synthetic_input(dimensions=(2, 1, 1), theta=0.0, emission=False):
    nx, ny, nz = dimensions
    with TemporaryDirectory() as directory, ExitStack() as stack:
        directory = Path(directory)
        input_path = directory / "plasma-vdsolverf.inp"
        thermal = 0.01 if emission else 1.0
        boundaries = "0,2,0" if emission else "2,0,0"
        emission_settings = (
            "nflag_emit=1, nepl=1, nemd=2, curf=1, curfs=15*-1, "
            "xmine=7.8, xmaxe=8, ymine=0.4, ymaxe=0.4, zmine=7.8, zmaxe=8"
            if emission else "nflag_emit=10*0, nepl=10*0"
        )
        spa, spe, speth = (0.1, np.sqrt(2.0), -np.pi / 4) if emission else (0, 0, 0)
        input_path.write_text(
            f"&esorem emflag=0 /\n"
            f"&plasma wp=10*0, wc=0, phiz={theta}, phixy={0 if emission else 31} /\n"
            f"&tmgrid nx={nx}, ny={ny}, nz={nz}, dt=0.1 /\n"
            f"&system nspec=1, npbnd={boundaries} /\n"
            f"&intp qm=10*1, peth={thermal}, path={thermal}, "
            f"spa={spa}, spe={spe}, speth={speth}, vdri=10*0, vdthz=10*0, vdthxy=10*0 /\n"
            f"&ptcond boundary_type='none' /\n"
            f"&emissn use_raycast=.false., {emission_settings} /\n"
        )
        data = SimpleNamespace(directory=directory, inp=SimpleNamespace(nx=nx, ny=ny, nz=nz))

        @contextmanager
        def temporary_input(data, suffix=None):
            yield SimpleNamespace(tmppath=input_path)

        # Stub file loading only: ctypes signatures and the Fortran solver are real.
        stack.enter_context(patch.object(wrapper.emout, "Emout", return_value=data))
        stack.enter_context(patch.object(wrapper, "TempolaryInput", temporary_input))
        stack.enter_context(patch.object(
            wrapper, "create_relocated_ebvalues",
            return_value=np.zeros((nz + 1, ny + 1, nx + 1, 9), dtype=np.float64),
        ))
        yield directory


@unittest.skipUnless(LIBRARY_PATH.is_file(), "requires a built Linux shared library")
class NativeRegressionTest(unittest.TestCase):
    def test_collision_time_at_nonunit_fixed_and_adaptive_steps(self):
        with synthetic_input() as directory:
            for adaptive in (False, True):
                for dt in (0.1, 0.2):
                    with self.subTest(adaptive=adaptive, dt=dt):
                        ts, probability, positions, velocities = wrapper.get_backtrace(
                            directory, 0, 0, Particle([1.95, 0.5, 0.5], [-2, 0, 0]),
                            dt, 3, use_adaptive_dt=adaptive, library_path=LIBRARY_PATH,
                        )
                        self.assertAlmostEqual(ts[-1], 0.025)
                        np.testing.assert_allclose(positions[-1], [2, 0.5, 0.5])
                        self.assertTrue(np.isfinite(probability))

    def test_oblique_isotropic_sources_remain_rotation_invariant(self):
        expected = np.exp(-0.5) / (2 * np.pi) ** 1.5
        for theta in (0.0, 30.0, 45.0, 90.0):
            with self.subTest(theta=theta), synthetic_input(theta=theta) as directory:
                ts, probability, positions, velocities = wrapper.get_backtrace(
                    directory, 0, 0, Particle([1, 0.5, 0.5], [-1, 0, 0]),
                    2.0, 3, library_path=LIBRARY_PATH,
                )
                self.assertAlmostEqual(probability, expected)

    def test_periodic_seam_deposits_only_visited_cells(self):
        with synthetic_input((8, 1, 8), emission=True) as directory:
            result = wrapper.estimate_velocity_range_map(
                directory, 0, 0, dt=0.2, max_step=1, coverage_sigma=1.0,
                safety_factor=1.0, minimum_count=1, collect_moments=True,
                show_progress=False, accumulator_cache_size=2, n_threads=1,
                library_path=LIBRARY_PATH,
            )
            visited = np.zeros_like(result.count, dtype=bool)
            visited[7, 0, 7] = True
            visited[0, 0, 0] = True
            self.assertGreater(result.count[7, 0, 7], 0)
            self.assertEqual(result.count[7, 0, 7], result.count[0, 0, 0])
            np.testing.assert_array_equal(result.count[~visited], 0)
            np.testing.assert_array_equal(result.weight_sum[~visited], 0)


if __name__ == "__main__":
    unittest.main()
