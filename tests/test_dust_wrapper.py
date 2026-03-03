import importlib
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch


class _DummyRotation:
    @staticmethod
    def from_euler(*args, **kwargs):
        class _Rot:
            def apply(self, vec):
                return vec

        return _Rot()


def _load_wrapper_module_once():
    if "vdsolverf.emses.wrapper" in sys.modules:
        return sys.modules["vdsolverf.emses.wrapper"]

    emout_stub = types.ModuleType("emout")
    emout_stub.Emout = object

    scipy_stub = types.ModuleType("scipy")
    spatial_stub = types.ModuleType("scipy.spatial")
    transform_stub = types.ModuleType("scipy.spatial.transform")
    transform_stub.Rotation = _DummyRotation
    spatial_stub.transform = transform_stub
    scipy_stub.spatial = spatial_stub

    with patch.dict(
        sys.modules,
        {
            "emout": emout_stub,
            "scipy": scipy_stub,
            "scipy.spatial": spatial_stub,
            "scipy.spatial.transform": transform_stub,
        },
    ):
        return importlib.import_module("vdsolverf.emses.wrapper")


WRAPPER = _load_wrapper_module_once()


class TestDustBacktraceWrapper(unittest.TestCase):
    def test_get_dust_backtrace_rejects_non_positive_max_step(self):
        with self.assertRaisesRegex(ValueError, "max_step"):
            WRAPPER.get_dust_backtrace_dll(
                directory=Path("."),
                istep=0,
                dust=object(),
                dt=1.0,
                max_step=0,
                use_adaptive_dt=False,
                max_probability_types=10,
                gravity=1.23,
                dll=object(),
            )

    def test_get_dust_backtrace_uses_legacy_os_alias(self):
        called = {}

        def fake_cdll(path):
            called["path"] = str(path)
            return object()

        def fake_impl(**kwargs):
            called["dll"] = kwargs["dll"]
            called["gravity"] = kwargs["gravity"]
            return "ok"

        with patch.object(WRAPPER, "CDLL", side_effect=fake_cdll), patch.object(
            WRAPPER, "get_dust_backtrace_dll", side_effect=fake_impl
        ):
            ret = WRAPPER.get_dust_backtrace(
                directory=Path("."),
                istep=0,
                dust=object(),
                dt=0.1,
                max_step=10,
                gravity=9.81,
                os="linux",
            )

        self.assertEqual(ret, "ok")
        self.assertTrue(called["path"].endswith("libvdist-solver-fortran.so"))
        self.assertEqual(called["gravity"], 9.81)


if __name__ == "__main__":
    unittest.main()
