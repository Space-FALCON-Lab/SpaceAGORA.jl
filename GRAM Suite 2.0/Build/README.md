# Using Make to build GRAMs

1. Build settings are auto-selected by host OS in `makefile.defs` (Darwin -> macOS profile, Linux -> gcc profile, Windows -> mingw profile).
2. make clean
3. make -j
4. Executables are in Build/bin.  Libraries are in Build/lib.
5. Edit and copy spice.txt into Build/bin before running the examples.

`Build/` is a generated output tree, not a portable artifact bundle. Do not commit or reuse `Build/bin`, `Build/lib`, `Build/mod`, or per-target object directories across macOS, Linux, and Windows; rebuild them natively on each host.

To build for a specific body (such as Neptune) issue the command "make Neptune -j".
To build a shared GRAM library for FFI workflows (Julia/Python/etc.), issue the command "make shared -j".

## macOS notes

The build files use GNU Make features. On macOS, use `gmake` (Homebrew `make`) instead of `/usr/bin/make`.

1. Ensure `SPICE_LIB` points to a macOS CSPICE static archive.
   - Defaults in `makefile.defs.macos`:
     - Apple Silicon: `common/cspice/lib/cspice_macos_arm64.a`
     - Intel Mac: `common/cspice/lib/cspice_macos_x86_64.a`
   - You can override at build time:
     - `gmake SPICE_LIB=/absolute/path/to/cspice.a`
2. Build:
   - Binaries/static libs: `gmake clean && gmake -j`
   - Shared library for Julia: `gmake clean && gmake shared -j`
3. Shared output:
   - macOS: `Build/lib/libGRAM.dylib`
   - Linux: `Build/lib/libGRAM.so`

## Portable setup (Intel Mac, Apple Silicon Mac, Linux, Windows MinGW)

Use the helper script to install/normalize CSPICE archive names used by the makefiles:

1. Run:
   - `./setup_cspice.sh`
2. Build:
   - macOS: `gmake clean && gmake shared -j`
   - Linux (GNU make): `make clean && make shared -j`
   - Windows MinGW (MSYS2): `mingw32-make clean && mingw32-make shared -j`

Notes:
- On macOS, `setup_cspice.sh` expects Homebrew `cspice` (`brew install cspice`).
- On Linux x86_64, the bundled `cspice_gcc85.a` is auto-normalized to `cspice_linux_x86_64.a`.
- On Linux arm64/aarch64, provide a native CSPICE archive or set `SPICE_LIB=/absolute/path/to/cspice.a`.
- On Windows MinGW, bundled `common/cspice/lib/cspice_mingw64.a` is used by default.
- If you switch operating systems in the same working tree, run `make clean` first or use the wrapper scripts in `simulation/GRAM/`, which record the host platform and force a clean rebuild when needed.

Files and Subfolders:
- googletest: Builds the googletest library used to run unit tests.
- GRAM: Builds an all inclusive GRAM library and example program.
- GRAMCommon: Builds the GRAM common framework library.
- Planet: Builds the library, examples, and PlanetGRAM program for the specified planet.
- makefile.defs: The configuration file for the make system (auto-selects host OS profile).
- makefile.defs.gcc: Sample configuration file for the GNU Compilers on Linux.
- makefile.defs.macos: Sample configuration file for clang compilers on macOS.
- makefile.defs.mingw: Sample configuration file for the MinGW Compilers on Windows.
- makefile.defs.intel: Sample configuration file for the Intel Compilers on Linux. (experimental)
- makefile.defs.clang: Sample configuration file for the clang Compilers on Windows. (experimental)
- makefile: The primary makefile for the make system.
- setup_cspice.sh: Host helper script to provision/normalize CSPICE archives.
- spice.txt: Namelist input specifying the SPICE data path.
- bin: Contains the programs created by the make system (folder generated during build process).
  + *Body*GRAM: The main analysis program for the specified Body model.
  + *Body*Examples: Quick examples using the specified Body model.
  + *Body*Example_Cpp: A C++ example for the specified Body model.
  + *Body*Example_C: A C example for the specified Body model.
  + *Body*Example_F: A Fortran example for the specified Body model.
- lib: Contains the libraries created by the make system (folder generated during build process).
  + lib*Body*: Library containing the specified Body model interface.
  + libGRAMCommon: Library containing the GRAM common framework interface.
  + libgtest: Library containing the GoogleTest unit testing interface.
  + libGRAM: Library combining all of the above libraries.
- mod: Contains the Fortran module files created by the make system (folder generated during build process).
