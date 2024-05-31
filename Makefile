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
        PLATFORM:=unknown
    endif
endif

BUILD_DIR=./build

all: $(LIBNAME)

$(LIBNAME): build copy_static shared copy_shared

.PHONY: build
build:
	fpm build --profile=release

copy_static:
	cp $(shell find ./build/*/ -type f -name lib$(LIBNAME).a) $(BUILD_DIR)

shared: shared_$(PLATFORM)

shared_linux:
	gfortran -shared -o $(BUILD_DIR)/lib$(LIBNAME).so -Wl,--whole-archive $(BUILD_DIR)/lib$(LIBNAME).a -Wl,--no-whole-archive

shared_darwin: 
	gfortran -dynamiclib -install_name ${BUILD_DIR}/lib$(LIBNAME).dylib -static-libgfortran -static-libquadmath -static-libgcc -o $(BUILD_DIR)/lib$(LIBNAME).dylib -Wl,-all_load $(BUILD_DIR)/lib$(LIBNAME).a -Wl,-noall_load

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
clean: 
	fpm clean --skip
	rm ${BUILD_DIR}/lib${LIBNAME}.a
	rm ${BUILD_DIR}/lib${LIBNAME}.so
