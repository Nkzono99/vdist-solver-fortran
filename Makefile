SHELL := /bin/bash

LIBNAME ?= vdist-solver-fortran
PYTHON ?= $(shell if [ -x .venv/bin/python ]; then echo $(CURDIR)/.venv/bin/python; \
	elif command -v python >/dev/null 2>&1; then echo python; \
	elif command -v python3.12 >/dev/null 2>&1; then echo python3.12; \
	elif command -v python3.11 >/dev/null 2>&1; then echo python3.11; \
	else echo python3; fi)
FPM ?= fpm
# Recovers fpm's incremental state after a failed or interrupted run.
FPM_GUARD ?= $(CURDIR)/scripts/fpm-stale-guard.sh
PROFILE ?= release
INSTALL_PROFILE ?= auto
PREFIX ?= $(PWD)
BATCH_FC ?= gfortran
MPI_FC ?= $(shell command -v mpif90 >/dev/null 2>&1 && echo mpif90 || echo mpiifx)
BATCH_FLAGS ?= $(if $(findstring ifx,$(BATCH_FC)),-O3 -qopenmp -fp-model precise,-O3 -fopenmp)
MPI_FLAGS ?= $(if $(or $(findstring mpiifx,$(MPI_FC)),$(findstring mpiifort,$(MPI_FC))),-O3 -qopenmp -fp-model precise,-O3 -fopenmp)
BATCH_PREFIX ?= $(CURDIR)/build/batch-serial
MPI_PREFIX ?= $(CURDIR)/build/batch-mpi

.PHONY: \
	all install install-local install-auto install-generic install-camphor \
	build test clean package-sdist package-check batch batch-mpi

all: install

install:
	BUILD_PROFILE=$(INSTALL_PROFILE) PREFIX=$(PREFIX) ./install.sh

install-local:
	BUILD_PROFILE=$(INSTALL_PROFILE) PREFIX=$(PWD)/.local ./install.sh

install-auto:
	BUILD_PROFILE=auto PREFIX=$(PREFIX) ./install.sh

install-generic:
	BUILD_PROFILE=generic PREFIX=$(PREFIX) ./install.sh

install-camphor:
	BUILD_PROFILE=camphor PREFIX=$(PREFIX) ./install.sh

build:
	$(FPM_GUARD) $(FPM) build --profile $(PROFILE) --flag "-fPIC -fopenmp"

test:
	$(FPM_GUARD) $(FPM) test --profile debug --flag "-fopenmp"

# Separate prefixes keep the serial executable available after an MPI build.
# These targets do not replace the Python shared library or change `all`.
batch:
	$(FPM_GUARD) $(FPM) install --profile $(PROFILE) --compiler $(BATCH_FC) --flag "$(BATCH_FLAGS)" --prefix "$(BATCH_PREFIX)"

batch-mpi:
	$(FPM_GUARD) $(FPM) install --profile $(PROFILE) --compiler $(MPI_FC) --flag "$(MPI_FLAGS) -DVDSOLVERF_MPI" --prefix "$(MPI_PREFIX)"

package-sdist:
	$(PYTHON) -m pip install --upgrade build
	rm -rf dist
	cd .. && $(PYTHON) -m build --sdist "$(CURDIR)" --outdir "$(CURDIR)/dist"

package-check: package-sdist
	$(PYTHON) -m pip install --upgrade twine
	$(PYTHON) -m twine check dist/*

clean:
	$(FPM) clean --skip
	rm -rf build dist *.egg-info lib
	rm -f vdsolverf/lib$(LIBNAME).so vdsolverf/lib$(LIBNAME).dylib vdsolverf/lib$(LIBNAME).dll
