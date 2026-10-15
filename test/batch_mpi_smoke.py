"""Reproducible, small MPI fixture; launch binaries separately in an allocation.

    .venv/bin/python test/batch_mpi_smoke.py prepare build/batch-smoke
    build/batch-serial/bin/vdsolverf-batch build/batch-smoke/serial/batch.nml
    mpiexec.hydra -n 4 build/batch-mpi/bin/vdsolverf-batch build/batch-smoke/mpi/batch.nml
    .venv/bin/python test/batch_mpi_smoke.py compare build/batch-smoke
"""
import json
from pathlib import Path
import sys

import f90nml
import numpy as np

# Running a file under test/ otherwise finds an older installed distribution.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vdsolverf.io import write_batch, read_batch_results


def prepare(root):
    bundle = root / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    nml = dict(
        esorem=dict(emflag=0), plasma=dict(wp=[0.]*10, wc=0., phiz=0., phixy=0.),
        tmgrid=dict(nx=3, ny=2, nz=4, dt=0.1), system=dict(nspec=1, npbnd=[[2, 2, 2]]*10),
        intp=dict(qm=[1.]*10, peth=[1.]*10, path=[1.]*10, vdri=[0.]*10, vdthz=[0.]*10,
                  vdthxy=[0.]*10, spa=[0.]*10, spe=[0.]*10, speth=[0.]*10),
        ptcond=dict(boundary_type="none"), emissn=dict(nepl=[0]*10, nflag_emit=[0]*10),
    )
    f90nml.Namelist(nml).write(bundle / "plasma-vdsolverf.inp", force=True)
    fields = np.arange(5*3*4*9, dtype="<f8").reshape(5, 3, 4, 9)*1e-5
    np.save(bundle / "ebvalues_ispec0.npy", fields)
    metadata = dict(schema="vdsolverf.field_bundle", version="0.1.0", nx=3, ny=2, nz=4,
                    directory=str(root.resolve()), istep=0, ispecs=[0], use_electric_field=True,
                    use_magnetic_field=True)
    (bundle / "bundle.json").write_text(json.dumps(metadata), encoding="utf-8")
    targets = np.array([[0.4+0.1*i, 0.8, 1.7, -4., 4., -4., 4., -4., 4.] for i in range(19)])
    targets[::5, 3:] = [-0.5, 0.5]*3
    targets[-1, 4] = -5.
    for name in ("serial", "mpi", "mpi-one"):
        write_batch(root / name, bundle, 0, targets, [-1., 0., 0.25, 1.],
                    dt=0.7, max_step=32, scout_bins=[3, 4, 3], max_depth=2,
                    max_samples_per_cell=4096, max_leaves_per_cell=128,
                    max_expansions=0, chunk_size=2, n_threads=4)
    print(f"Prepared 19 targets / 10 chunks in {root}")


def compare(root):
    serial = read_batch_results(root / "serial")
    mpi = read_batch_results(root / "mpi")
    assert mpi.metadata["summary"]["ranks"] == 4, mpi.metadata["summary"]
    for name in ("positions", "psi_grid", "flux", "flux_error", "status", "sample_count", "terminal_box_count"):
        np.testing.assert_array_equal(getattr(mpi, name), getattr(serial, name))
    assert serial.metadata["summary"]["status_counts"] == mpi.metadata["summary"]["status_counts"]
    print("PASS: 4-rank MPI and serial outputs are bit-identical (all arrays and status counts)")
    print(json.dumps(mpi.metadata["summary"], indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("prepare", "compare"):
        raise SystemExit("usage: batch_mpi_smoke.py prepare|compare ROOT")
    globals()[sys.argv[1]](Path(sys.argv[2]).resolve())
