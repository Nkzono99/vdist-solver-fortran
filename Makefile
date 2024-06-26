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
export FPM_FFLAGS := -Ofast -m64 -fPIC -fopenmp -lgomp
export FPM_FFLAGS := $(FPM_FFLAGS) -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow

BUILD_DIR=./build

all: $(LIBNAME)

$(LIBNAME): build copy_static shared copy_shared

.PHONY: build
build:
	fpm build --profile=release --verbose

copy_static:
	cp $(shell find ./build/*/ -type f -name lib$(LIBNAME).a) $(BUILD_DIR)

shared: shared_$(PLATFORM)

shared_linux:
	gfortran -shared -o $(BUILD_DIR)/lib$(LIBNAME).so -Wl,--whole-archive $(BUILD_DIR)/lib$(LIBNAME).a -Wl,--no-whole-archive

shared_darwin: 
	gfortran ${BUILD_DIR}/lib$(LIBNAME).a -dynamiclib -install_name ${BUILD_DIR}/lib$(LIBNAME).dylib -static-libgfortran -static-libquadmath -static-libgcc -o $(BUILD_DIR)/lib$(LIBNAME).dylib -Wl,-all_load $(BUILD_DIR)/lib$(LIBNAME).a -Wl,-noall_load

shared_windows: 
	gfortran -shared -static -o $(BUILD_DIR)/lib$(LIBNAME).dll -Wl,--out-implib=$(BUILD_DIR)/lib$(LIBNAME).dll.a,--export-all-symbols,--enable-auto-import,--whole-archive $(BUILD_DIR)/lib$(LIBNAME).a -Wl,--no-whole-archive

copy_shared: copy_shared_${PLATFORM}

copy_shared_linux:
	cp ${BUILD_DIR}/lib${LIBNAME}.so vdsolverf/

copy_shared_darwin:
	cp ${BUILD_DIR}/lib${LIBNAME}.dylib vdsolverf/

copy_shared_windows:
	cp ${BUILD_DIR}/lib${LIBNAME}.dll vdsolverf/

.PHONY: clean
clean: clean_$(PLATFORM)
	fpm clean --skip
	rm ${BUILD_DIR}/lib${LIBNAME}.a

clean_linux:
	rm vdsolverf/lib${LIBNAME}.so

clean_darwin:
	rm vdsolverf/lib${LIBNAME}.dylib

clean_windows:
	rm vdsolverf/lib${LIBNAME}.dll
