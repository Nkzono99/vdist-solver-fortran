"""Reject EMSES source settings that the static vdist replay cannot reproduce.

vdist rebuilds the EMSES boundary and surface sources from the namelist and a
field snapshot. Settings that EMSES applies but vdist would ignore are rejected
here, from the original EMSES namelist, before the temporary vdist input drops
them. The native ``validate_emses_input`` checks the remaining rules on the
temporary input. See ``docs/emses_sync.md``.
"""
import numpy as np


_UNSUPPORTED_NONZERO = {
    "periodic_raycast": "periodic raycast shadows",
    "omniemit": "omnidirectional patch emission",
    "boundary_mirror_reflection_rate": "surface mirror reflection",
    "boundary_reversal_reflection_rate": "surface reversal reflection",
    "boundary_cosine_reflection_rate": "surface cosine reflection",
    "enable_secondary_electron_emission": "secondary electron emission",
}


def validate_source_support(namelist):
    """Raise ``ValueError`` for source physics outside the static replay.

    ``namelist`` is the parsed EMSES namelist (groups of key/value pairs).
    """
    values = {key.lower(): value for group in namelist.values() if hasattr(group, "items")
              for key, value in group.items()}

    def array(key, default=0):
        def fill_missing(value):
            if isinstance(value, (list, tuple)):
                return [fill_missing(item) for item in value]
            return default if value is None else value
        return np.asarray(fill_missing(values.get(key, default)), dtype=float)

    for key, name in _UNSUPPORTED_NONZERO.items():
        if np.any(array(key) != 0):
            raise ValueError(f"{key}: {name} is not supported by the static EMSES source replay")

    if np.any(array("f", 1) != 1):
        raise ValueError("f: loss-cone/ring velocity distributions are not supported by the static EMSES source replay")

    # EMSES splits surface-emission tables into two thermal components when an
    # entry beyond nspec is positive and mainprop < 1.
    nspec = int(values.get("nspec", 0))
    if nspec > 0 and float(values.get("mainprop", 1.0)) < 1.0:
        for key in ("path", "peth"):
            if np.any(array(key).ravel()[nspec:] > 0):
                raise ValueError(f"{key}(nspec+1): two-component surface emission (mainprop < 1) "
                                 "is not supported by the static EMSES source replay")
