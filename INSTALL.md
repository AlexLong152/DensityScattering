# Installing the build dependencies

The kernel Fortran code links against parallel HDF5 (with Fortran
bindings) and reads density files compressed with the
[h5z-zfp](https://github.com/LLNL/H5Z-ZFP) plugin. Both must be
installed before the kernels will build and read densities.

If your system package manager has up-to-date versions
(e.g. `apt install libhdf5-openmpi-dev` on recent Debian/Ubuntu),
that is the fast path. The rest of this guide covers source builds
for the common case where the packaged version is too old or missing
the bindings the project needs.

## 1. HDF5 with OpenMPI

Instructions below are for Debian-based systems (Ubuntu, Linux Mint).
For Arch-based distributions see the note at the end of this section.

### 1.1 Install the required system packages

```bash
sudo apt install build-essential openmpi-bin libopenmpi-dev
```

### 1.2 Build and install HDF5 with OpenMPI

Pick a release from the [HDF5 download
page](https://support.hdfgroup.org/downloads/hdf5/), e.g. 1.14.6:

```bash
wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-<version>.tar.gz
tar -xvzf hdf5-<version>.tar.gz
cd hdf5-<version>
mkdir build && cd build

CC=mpicc FC=mpif90 ../configure \
  --enable-parallel \
  --enable-fortran \
  --prefix=/usr/local/hdf5/openmpi

make -j$(nproc)
sudo make install
```

This installs HDF5 system-wide under `/usr/local/hdf5/openmpi`, which
is the default `HDF5_DIR` in the repo's `config.mk`. If you install
elsewhere, edit `config.mk` (or set `HDF5_DIR` in your shell).

### 1.3 Runtime configuration

Add to `~/.bashrc` (or your shell's equivalent):

```bash
export LD_LIBRARY_PATH=/usr/local/hdf5/openmpi/lib:$LD_LIBRARY_PATH
export PATH=/usr/local/hdf5/openmpi/bin:$PATH
```

The kernel Makefiles also embed an `RPATH` at link time (via
`-Wl,-rpath,$(HDF5_DIR)/lib`), so the executables find the HDF5
shared libraries without `LD_LIBRARY_PATH` once they're built — but
the `LD_LIBRARY_PATH` line is needed for any other tool that links
against HDF5.

### Note for Arch-based distributions

On Manjaro / Arch, `pacman -S hdf5-openmpi` provides HDF5 itself, but
conflicts with `pacman -S hdf5` if both are selected. If the spack
install of h5z-zfp in section 2 fails to find a usable HDF5, build
HDF5 from source as above.

## 2. h5z-zfp compression plugin
This method is NOT optimal, and might not work on all systems. Follow these instructions at your own risk.
I believe this was required for Aarch based systems, but Debian based systems might just work with the `hdf5-plugin-lzp` or similar package, which is much easier to install. What is below should function though at least as a fallback.

The density files distributed by the Jülich store are compressed with
[ZFP](https://github.com/LLNL/zfp) and require the `h5z-zfp` HDF5
filter plugin to read. Without it, any kernel run will abort with:

```
H5PLpath.c line 857 in H5PL__find_plugin_in_path(): can't open directory
(/usr/local/hdf5/lib/plugin). Please verify its existence
```

The recommended install path is via [spack](https://spack.io/).

### 2.1 Install spack

```bash
git clone https://github.com/llnl/spack.git
. spack/share/spack/setup-env.sh
```

### 2.2 Install h5z-zfp

```bash
spack install h5z-zfp
```

This compile takes several minutes.

### 2.3 Copy the plugin into HDF5's plugin directory

Find where spack put the plugin:

```bash
spack find -vp h5z-zfp
```

The output looks like
`<spack-root>/opt/spack/<arch>/<compiler>/h5z-zfp-develop-<hash>`. That
directory contains a `plugin/` subdirectory; copy it into HDF5's
plugin directory (`/usr/local/hdf5/lib/plugin/` for the default
install above):

```bash
sudo cp -r <h5z-zfp-prefix>/plugin /usr/local/hdf5/lib/
```

If `/usr/local/hdf5/lib/plugin/` already contains other HDF5 plugins,
copy `h5z-zfp` into a subdirectory (e.g.
`/usr/local/hdf5/lib/plugin/h5z-zfp/`) to avoid overwriting them.
