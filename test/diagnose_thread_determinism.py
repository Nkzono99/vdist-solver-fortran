"""Raw C-ABI reproducibility probe; run on a compute node, not unittest discovery.

run LIBRARY OUTPUT_DIRECTORY evaluates the original nonzero-field batch fixture
three times using OMP_NUM_THREADS. compare OUTPUT_DIRECTORY compares raw bits
against t1-r0, reporting every changed output, spatial point and max abs/rel error.
"""
import ctypes as ct
import os
from pathlib import Path
import sys

import numpy as np


def run(library, root):
    root.mkdir(parents=True, exist_ok=True)
    # Matches test/batch_mpi_smoke.py, including an invalid final target.
    path = root / "plasma.inp"
    path.write_text("""&esorem emflag=0 /
&plasma wp=10*0, wc=0, phiz=0, phixy=0 /
&tmgrid nx=3, ny=2, nz=4, dt=0.1 /
&system nspec=1, npbnd=30*2 /
&intp qm=10*1, peth=10*1, path=10*1, vdri=10*0, vdthz=10*0,
 vdthxy=10*0, spa=10*0, spe=10*0, speth=10*0 /
&ptcond boundary_type='none' /
&emissn nepl=10*0, nflag_emit=10*0 /
""")
    fields = np.arange(5*3*4*9, dtype=np.float64).reshape(5, 3, 4, 9)*1e-5
    points = np.array([[0.4+0.1*i, 0.8, 1.7] for i in range(19)])
    bounds = np.tile([-4., 4., -4., 4., -4., 4.], (19, 1))
    bounds[::5] = [-0.5, 0.5]*3
    bounds[-1, 1] = -5.
    psi = np.array([-1., 0., 0.25, 1.])
    grid = np.array(np.meshgrid(np.linspace(-4, 4, 3), np.linspace(-4, 4, 4),
                               np.linspace(-4, 4, 3), indexing="ij")).reshape(3, -1).T
    velocities = np.ascontiguousarray(np.tile(grid, (19, 1)))
    positions = np.ascontiguousarray(np.repeat(points, len(grid), axis=0))
    dll = ct.CDLL(str(library.resolve()))
    I, D, P = ct.c_int, ct.c_double, ct.c_void_p
    common_types = [ct.c_char_p, I, I, I, I, P, I]
    dll.evaluate_dust_currents.argtypes = common_types + [I, P, P, I, P, D] + [I]*9 + [D]*3 + [I] + [P]*6
    dll.get_probabilities_octree_with_terminals.argtypes = common_types + [I, P, P, D] + [I]*9 + [D]*3 + [I] + [P]*17
    dll.get_probabilities.argtypes = common_types + [I, P, P, D, I, I, I] + [P]*4
    dll.estimate_velocity_range_map.argtypes = common_types + [D]*3 + [I]*9 + [P]*13
    for name in ("evaluate_dust_currents", "get_probabilities_octree_with_terminals",
                 "get_probabilities", "estimate_velocity_range_map"):
        getattr(dll, name).restype = None
    def ptr(a):
        return a.ctypes.data_as(P)
    threads = I(int(os.environ.get("OMP_NUM_THREADS", "1")))
    encoded = str(path.resolve()).encode()
    common = [encoded, len(encoded), 3, 2, 4, ptr(fields), 1]
    for repeat in range(3):
        out = {}
        def array(name, shape, integer=False):
            out[name] = np.zeros(shape, dtype=np.int32 if integer else np.float64)
            return ptr(out[name])
        dll.evaluate_dust_currents(*common, 19, ptr(points), ptr(bounds), 4, ptr(psi), 0.7,
            32, 0, 100, 3, 4, 3, 2, 4096, 128, 1e-4, 1e-4, 1.5, 0,
            array("dust_flux", (19, 4)), array("dust_error", (19, 4)),
            array("dust_status", 19, True), array("dust_samples", 19, True),
            array("dust_terminals", 19, True), ct.byref(threads))
        dll.get_probabilities_octree_with_terminals(*common, 19, ptr(points), ptr(bounds), 0.7,
            32, 0, 100, 3, 4, 3, 2, 4096, 128, 1e-4, 1e-4, 1.5, 0,
            array("oct_sample_spatial", (19, 4096), True), array("oct_velocity", (19, 4096, 3)),
            array("oct_probability", (19, 4096)), array("oct_leaf_spatial", (19, 128), True),
            array("oct_bounds", (19, 128, 6)), array("oct_min", (19, 128)), array("oct_max", (19, 128)),
            array("oct_depth", (19, 128), True), array("oct_start", (19, 128), True),
            array("oct_count", (19, 128), True), array("oct_status", 19, True),
            array("oct_samples", 19, True), array("oct_leaves", 19, True),
            array("oct_total_samples", 1, True), array("oct_total_leaves", 1, True),
            ct.byref(threads), array("oct_terminal", (19, 128), True))
        n = len(positions)
        dll.get_probabilities(*common, n, ptr(positions), ptr(velocities), 0.7, 32, 0, 100,
            array("prob_probability", (19, len(grid))), array("prob_position", (19, len(grid), 3)),
            array("prob_velocity", (19, len(grid), 3)), ct.byref(threads))
        dll.estimate_velocity_range_map(*common, 0.7, 3., 1.1, 16, 0, 100, 1, 0, 1, 1, 0, 4,
            *[array("range_"+name, (4, 2, 3)) for name in ("vx_min", "vx_max", "vy_min", "vy_max", "vz_min", "vz_max")],
            array("range_count", (4, 2, 3), True), array("range_weight", (4, 2, 3)),
            array("range_mean", (4, 2, 3, 3)), array("range_cov", (4, 2, 3, 3, 3)),
            array("range_status", (4, 2, 3), True), array("range_confidence", (4, 2, 3)), ct.byref(threads))
        target = root / f"t{threads.value}-r{repeat}.npz"
        np.savez(target, **out)
        print("saved", target, flush=True)


def compare(root):
    serial = root / "t1-r0.npz"
    for path in sorted(root.glob("t*-r*.npz")):
        values = np.load(path)
        compare_pair(serial, path, np.load(serial), values)
        same_threads = root / (path.stem.split("-")[0] + "-r0.npz")
        if same_threads != serial and same_threads != path:
            compare_pair(same_threads, path, np.load(same_threads), values)


def compare_pair(reference, path, baseline, values):
    label = f"{path.name} vs {reference.name}"
    changed = []
    for key in baseline.files:
        a, b = baseline[key], values[key]
        bits_a = a.view(np.uint64 if a.dtype.itemsize == 8 else np.uint32)
        bits_b = b.view(bits_a.dtype)
        mask = bits_a != bits_b
        if not mask.any():
            continue
        changed.append(key)
        delta = np.abs(a.astype(float) - b.astype(float))
        relative = delta / np.maximum(np.abs(a), np.finfo(float).tiny)
        # Range coordinates are grid (z,y,x); other leading indices are spatial points.
        indices = np.argwhere(mask)
        print(label, key, "unequal", int(mask.sum()), "max_abs", float(np.nanmax(delta)),
              "max_rel", float(np.nanmax(relative)), "indices", indices[:30].tolist(),
              "points", np.unique(indices[:, 0]).tolist())
    print(label, "changed:", changed or "NONE", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(Path(sys.argv[2]), Path(sys.argv[3]))
    elif sys.argv[1] == "compare":
        compare(Path(sys.argv[2]))
    else:
        raise SystemExit("usage: diagnose_thread_determinism.py run LIBRARY ROOT | compare ROOT")
