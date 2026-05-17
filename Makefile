SHELL := /bin/bash

LIBNAME ?= vdist-solver-fortran
PYTHON ?= $(shell if [ -x .venv/bin/python ]; then echo $(CURDIR)/.venv/bin/python; \
	elif command -v python >/dev/null 2>&1; then echo python; \
	elif command -v python3.12 >/dev/null 2>&1; then echo python3.12; \
	elif command -v python3.11 >/dev/null 2>&1; then echo python3.11; \
	else echo python3; fi)
FPM ?= fpm
PROFILE ?= release
INSTALL_PROFILE ?= auto
PREFIX ?= $(PWD)

.PHONY: \
	all install install-local install-auto install-generic install-camphor \
	build test clean package-sdist package-check

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
	$(FPM) build --profile $(PROFILE) --flag "-fPIC -fopenmp"

test:
	$(FPM) test --profile debug --flag "-fopenmp"

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
