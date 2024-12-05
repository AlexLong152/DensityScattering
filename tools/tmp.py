import os

import pandas as pd
from nucdens import access
import numpy as np
from pyperclip import copy

pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)

# connect to database, choose working directory for downloading densities
# workdir = os.environ["HOME"]+"/work/densitywork"
workdir = os.environ["HOME"] + r"/EFFNUCLEON/COMPTON/FEW-NUCLEON/DENSITIES/"
workdir = os.environ["HOME"] + r"/OneDrive/"

beta_webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore-beta/"
webbase = "https://just-object.fz-juelich.de:8080/v1/AUTH_1715b19bd3304fb4bb04f4beccea0cf2/densitystore/"
orig_webbase = "https://datapub.fz-juelich.de/anogga/files/densindx/"
try:
    densdf = access.database(workdir=workdir, webbase=beta_webbase)
    # densdf = access.database(workdir=workdir, webbase=webbase)
except BaseException:
    # use mirror if server is offline
    print("Using mirror database")
    densdf = access.database(workdir=workdir, webbase=orig_webbase)

print("Proceeding with workdir=" + workdir)
