# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import runfolderTwobody as rt


# 3He
# densfolder = r"/home/alexander/OneDrive/densities-3He/2Ndensities/132MeV/60deg/"
# outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/2bod/132MeV/"

# 6Li
# densfolder = r"/home/alexander/OneDrive/densities-6Li/2Ndensities/133MeV/60deg/"
# outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/2bod/"

# 4He
# densfolder = r"/home/alexander/OneDrive/densities-4He/2Ndensities/133MeV/60deg"
# outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/2bod/"

# 3H
densfolder = "/home/alexander/OneDrive/densities-3H/2Ndensities/132MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/2bod"
rt.Odelta = 4
rt.j12max = 1
rt.batch = 1
rt.executableName = r"run.twobodyvia2Ndensity.PionPhotoProdThresh"
rt.main(densfolder, outputfolder)
