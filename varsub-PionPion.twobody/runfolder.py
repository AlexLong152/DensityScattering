# -*- coding: utf-8 -*-

"""
@author: alexl
"""

# import numpy as np
import runfolderTwobody as rt


densfolder = r"/home/alexander/OneDrive/densities-3He/2Ndensities/132MeV/0deg/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-3He/2bod/thresh/"
rt.Odelta = 4
rt.j12max = 2

rt.executableName = r"run.twobodyvia2Ndensity.PionPion"
rt.main(densfolder, outputfolder)
