# -*- coding: utf-8 -*-

"""
@author: alexl
"""

# import numpy as np
import runfolderTwobody as rt


densfolder = r"/home/alexander/OneDrive/densities-3He/2Ndensities/136MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-3He/2bod/chargedThresh/"
rt.Odelta = 4
rt.j12max = 1

rt.executableName = r"run.twobodyvia2Ndensity.PionPion"
rt.main(densfolder, outputfolder)
