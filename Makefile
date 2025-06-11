SHELL=/bin/bash

LIBNAME=vdist-solver-fortran

ifeq ($(OS),Windows_NT)
    PLATFORM := windows
else
    UNAME_OS := $(shell uname -s)
    ifeq ($(UNAME_OS),Linux)
        PLATFORM := linux
    else ifeq ($(UNAME_OS),Darwin)
        PLATFORM := darwin
    else
        PLATFORM := unknown
    endif
endif

export FPM_FC := gfortran
export FPM_FFLAGS := -Ofast -m64 -fPIC -fopenmp
export FPM_FFLAGS := $(FPM_FFLAGS) -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
BUILD_DIR=build

all: $(LIBNAME)

$(LIBNAME): build shared copy_shared

.PHONY: build
build: build_$(PLATFORM)

build_linux:
	fpm install --profile=release --prefix ./

build_darwin:
	fpm install --profile=release --archiver /usr/bin/ar --prefix ./

build_windows:
	fpm install --profile=release --prefix ./

shared: shared_$(PLATFORM)

shared_linux:
	gfortran -shared -o lib/lib$(LIBNAME).so -Wl,--whole-archive lib/lib$(LIBNAME).a -Wl,--no-whole-archive -fopenmp

shared_darwin: 
	gfortran lib/lib$(LIBNAME).a -dynamiclib -install_name lib/lib$(LIBNAME).dylib -static-libgfortran -static-libquadmath -static-libgcc -o lib/lib$(LIBNAME).dylib -Wl,-all_load lib/lib$(LIBNAME).a -Wl,-noall_load -fopenmp

shared_windows: 
	gfortran -shared -static -o lib\\lib$(LIBNAME).dll -Wl,--out-implib=lib\\lib$(LIBNAME).dll,--export-all-symbols,--enable-auto-import,--whole-archive lib/lib$(LIBNAME).a -Wl,--no-whole-archive -fopenmp

copy_shared: copy_shared_${PLATFORM}

copy_shared_linux:
	cp lib/lib${LIBNAME}.so vdsolverf/

copy_shared_darwin:
	cp lib/lib${LIBNAME}.dylib vdsolverf/

copy_shared_windows:
	copy /y lib\\lib${LIBNAME}.dll vdsolverf\\lib${LIBNAME}.dll

.PHONY: clean
clean: clean_$(PLATFORM)

clean_linux:
	fpm clean --skip
	rm lib/lib${LIBNAME}.a
	rm vdsolverf/lib${LIBNAME}.so

clean_darwin:
	fpm clean --skip
	rm lib/lib${LIBNAME}.a
	rm vdsolverf/lib${LIBNAME}.dylib

clean_windows:
	fpm clean --skip
	del /Q lib/lib${LIBNAME}.a >NUL 2>NUL || echo ok
	del /Q "vdsolverf/lib${LIBNAME}.dll" >NUL 2>NUL || echo ok
