# config.mk -- shared build configuration.
#
# Every kernel Makefile (PionPion.onebody/, PionPhotoProd.onebody/,
# PionPhotoProdThresh.onebody/, varsub-PionPion.twobody/,
# varsub-PionPhotoProdThresh.twobody/) does `include ../config.mk`,
# so editing this file once propagates to all kernel builds.
#
# Path to your HDF5 installation (the directory containing include/ and lib/).
# `?=` means: if HDF5_DIR is already set in the environment, that wins.
# To find your HDF5 prefix on a typical system:
#   pkg-config --variable=prefix hdf5-openmpi-fortran

HDF5_DIR ?= /usr/local/hdf5/openmpi
