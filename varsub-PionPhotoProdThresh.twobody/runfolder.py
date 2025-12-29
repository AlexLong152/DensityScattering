# -*- coding: utf-8 -*-

"""
@author: alexl
"""

# import numpy as np
import runfolderTwobody as rt


# 4He
densfolder = r"/home/alexander/OneDrive/densities-4He/2Ndensities/133MeV/60deg/"
outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/"
rt.Odelta = 4
rt.j12max = 1
rt.batch = 1
rt.executableName = r"run.twobodyvia2Ndensity.PionPhotoProdThresh"
rt.main(densfolder, outputfolder)
