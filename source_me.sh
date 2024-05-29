export FPM_FFLAGS="-Ofast -m64 -fPIC"
# export FPM_FFLAGS="-fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow"
export FPM_FFLAGS="$FPM_FFLAGS -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow"
