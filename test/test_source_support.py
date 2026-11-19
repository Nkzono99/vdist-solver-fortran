import unittest

from vdsolverf.emses.source_support import validate_source_support


class SourceSupportTest(unittest.TestCase):
    def test_standard_sources_are_accepted(self):
        validate_source_support({
            "system": {"nspec": 2},
            "intp": {"path": [1.0, 2.0, 0.5], "peth": [1.0, 2.0], "f": [None, 1]},
            "emissn": {"pe_ray_cast": True, "periodic_raycast": False, "mainprop": 1.0},
        })

    def test_unsupported_physics_is_not_silently_discarded(self):
        for group in ({"periodic_raycast": True}, {"omniemit": [0, 1]}, {"f": [0.5]},
                      {"boundary_mirror_reflection_rate": [[0.2]]},
                      {"boundary_cosine_reflection_rate": [[0.1]]},
                      {"enable_secondary_electron_emission": [[True]]},
                      {"path": [1.0, 2.0], "mainprop": 0.8}):
            with self.subTest(group=group), self.assertRaises(ValueError):
                validate_source_support({"system": {"nspec": 1}, "parameters": group})

    def test_extra_thermal_entries_only_matter_with_a_mixture(self):
        # EMSES uses path(nspec+1) only when mainprop < 1.
        validate_source_support({"system": {"nspec": 1}, "intp": {"path": [1.0, 2.0]}})


if __name__ == "__main__":
    unittest.main()
