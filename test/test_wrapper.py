import tempfile
import unittest
from pathlib import Path
from typing import Tuple, get_type_hints
from unittest.mock import patch

import numpy as np

from vdsolverf.core import Particle
from vdsolverf.emses import wrapper


class _FakeInput:
    nx = 1
    ny = 1
    nz = 1


class _FakeData:
    inp = _FakeInput()


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


class _FakeDll:
    def __init__(self):
        self.get_backtraces = _FakeBacktracesFunction()


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


if __name__ == "__main__":
    unittest.main()
