# -*- coding: utf-8 -*-

"""
@author: alexl
"""

# import numpy as np
import runfolderTwobody as rt

# 3He settings
densfolder = r"/home/alexander/OneDrive/densities-3He/2Ndensities/136MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-3He/2bod/thresh/"

# 4He settings
outputfolder = r"/home/alexander/Dropbox/PionPion/results-4He/2bod/thresh/"
densfolder = r"/home/alexander/OneDrive/densities-4He/2Ndensities/136MeV/"

# 3H settings
densfolder = r"/home/alexander/OneDrive/densities-3H/2Ndensities/138MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-3H/2bod/thresh/"

# 6Li settings
densfolder = r"/home/alexander/OneDrive/densities-6Li/2Ndensities/138MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-6Li/2bod/thresh/"

rt.batch = 3
rt.Odelta = 4
rt.j12max = 1

rt.executableName = r"run.twobodyvia2Ndensity.PionPion"
rt.main(densfolder, outputfolder)
