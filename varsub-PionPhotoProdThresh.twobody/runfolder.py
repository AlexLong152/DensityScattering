# -*- coding: utf-8 -*-

"""
@author: alexl
"""

# import numpy as np
import runfolderTwobody as rt


densfolder = r"/home/alexander/OneDrive/densities-3He/2Ndensities/132MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/"
rt.Odelta = 4
rt.j12max = 2

rt.executableName = r"run.twobodyvia2Ndensity.PionPhotoProdThresh"
rt.main(densfolder, outputfolder)
