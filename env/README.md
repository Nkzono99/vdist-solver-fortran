# Build Profiles

`install.sh` loads host-specific build presets from this directory.

- `common.env`: shared defaults
- `camphor.env`: optimized Intel/OpenMP flags for camphor
- `profile.template.env`: template for new host profiles

Profile selection:

- `BUILD_PROFILE=auto` (default): detect by hostname (`camphor* -> camphor`, otherwise `generic`)
- `BUILD_PROFILE=generic`: portable gfortran/OpenMP build
- `BUILD_PROFILE=<name>`: load `env/<name>.env`

Variables in each `*.env` file use `: "${VAR:=...}"` style so user-provided
environment variables can still override them.

To add a new host profile:

1. `cp env/profile.template.env env/<profile>.env`
2. Fill `FC`, `FFLAGS`, `SHARED_LDFLAGS`, and module settings
3. Run `BUILD_PROFILE=<profile> ./install.sh` or
   `make install INSTALL_PROFILE=<profile>`
