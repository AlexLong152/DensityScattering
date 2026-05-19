import os
import sys
import warnings
from pandas.errors import PerformanceWarning
import numpy as np
import pandas as pd
from nucdens import access

FILENAME = "3He-onebody-piphoto-example.h5"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
for stale in ("example.dat", "densities_table.h5", "out.dat"):
    stale_path = os.path.join(SCRIPT_DIR, stale)
    if os.path.exists(stale_path):
        os.remove(stale_path)

REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
BUILD_DIRS = [
    os.path.join(REPO_ROOT, "common-densities", "density-modules"),
    os.path.join(REPO_ROOT, "common-densities", "varsub-density-modules"),
    os.path.join(REPO_ROOT, "PionPhotoProdThresh.onebody"),
]
RUN_DIR = os.path.join(REPO_ROOT, "PionPhotoProd.onebody")

INPUT_DAT = (
    "131.85 131.85 10                           omegaLow, omegaHigh, omegaStep\n"
    "60 60 15                             thetaLow, thetaHigh, thetaStep\n"
    "../example/out.dat\n"
    "'../example/3He-onebody-piphoto-example.h5'\n"
    "cm_symmetry_verbos_extQnumlimit=3                    frame, symmetry, verbosity of STDOUT\n"
    "Odelta2                                  Calctype -- VaryAXN varies NucelonN's amplitude AX (N=p,n; X=1-6)\n"
    "32\t\t\t             number of points Nx in Feynman parameter integration\n"
    "COMMENTS:_v2.0\n"
)


def error_compile(path):
    print("=" * 70)
    print(f"ERROR: `make` failed in {path}")
    print(
        "Please edit the Makefile (and make.inc, where applicable) for your\n"
        "system's compilers and library paths, then consult README.md and\n"
        "INSTALL.md for build instructions."
    )
    print("=" * 70)
    sys.exit(1)


def build(path):
    print(f"\n>>> Building in {path}")
    cwd = os.getcwd()
    try:
        os.chdir(path)
        rc = os.system("make")
    finally:
        os.chdir(cwd)
    if rc != 0:
        error_compile(path)


for d in BUILD_DIRS:
    build(d)


warnings.simplefilter(action="ignore", category=PerformanceWarning)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


# open database
beta_webbase = "https://just-object.fz-juelich.de:9000/jikp03/densitystore-beta/"
workdir = SCRIPT_DIR
target = os.path.join(workdir, FILENAME)
if os.path.exists(target):
    os.remove(target)

densdb = access.database(workdir=workdir, webbase=beta_webbase)
df = densdb.pddf

eps = 2

thetas = np.array([60])
theta_mask = [(df["theta"] - theta).abs() < eps for theta in thetas]
theta_mask = np.logical_or.reduce(theta_mask)

energies = np.array([133])
omega_masks = [(df["omega"] - e).abs() < eps for e in energies]
omega_mask = np.logical_or.reduce(omega_masks)

select = (
    (df["N"] == 1)
    & (df["Z"] == 2)
    & (df["kind"] == "one")
    & (df["LambdaNN"] == 450)
    & theta_mask
    & omega_mask
)
df = df[select]
if df.empty:
    raise RuntimeError("No densities matched the selection.")

row = df.iloc[0]
hashname, _ = densdb.get_file(hashname=row["hashname"], kind=row["kind"])
assert os.path.isfile(hashname)
os.replace(hashname, target)
print(f"Downloaded file to: {target}")

input_dat_path = os.path.join(RUN_DIR, "input.dat")
with open(input_dat_path, "w") as f:
    f.write(INPUT_DAT)
print(f"Wrote input file: {input_dat_path}")

print(f"\n>>> Running ./run.sh in {RUN_DIR}")
cwd = os.getcwd()
try:
    os.chdir(RUN_DIR)
    rc = os.system("./run.sh")
finally:
    os.chdir(cwd)
if rc != 0:
    print(f"run.sh exited with non-zero status ({rc}).")
    sys.exit(1)
