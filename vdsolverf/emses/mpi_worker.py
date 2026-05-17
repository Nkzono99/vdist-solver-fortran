"""Command-line worker for :mod:`vdsolverf.emses.mpi` launchers."""

from .mpi import worker_main


if __name__ == "__main__":
    raise SystemExit(worker_main())
