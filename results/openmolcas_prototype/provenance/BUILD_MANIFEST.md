# Build and host manifest

## Source

- Repository: <https://github.com/OpenMolcas/OpenMolcas>
- Release tag: `v26.06`
- Commit: `8355057f32d65706a35996b5ab07cac2962bb728`
- Submodules: initialized recursively; bundled LAPACK v3.12.0 at commit `04b044e020a3560ccfa9988c8a80a1fb7083fc2e`

## Toolchain

- CMake 3.31.10
- Ninja 1.11.1 (Kitware jobserver build)
- GNU Fortran 13.3.0 (`gfortran-13`, Ubuntu 24.04 package)
- OpenMolcas pymolcas driver py2.32

## Selected CMake cache values

```text
BUILD_TESTING=OFF
CMAKE_BUILD_TYPE=Release
GA=OFF
HDF5=OFF
LINALG=Internal
MPI=OFF
TOOLS=OFF
```

The source tree was clean at the recorded commit. The build completed to 100%, including `rasscf.exe`, `caspt2.exe`, `rassi.exe`, and `single_aniso.exe`; a no-change rebuild check returned all targets built.

## Host allocation

```text
Kernel: Linux 6.18.35 x86_64
CPU: Intel Xeon Platinum 8370C @ 2.80 GHz
Allocated CPUs: 9 (tests used 1 process and 1 thread)
RAM: 21 GiB
Swap: none
Available disk before runs: approximately 49 GiB
```

Execution date: 16 August 2026, America/New_York.
