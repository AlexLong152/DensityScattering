# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import runfolderOnebody as ro
from os import system

system("make clean;make")
ro.Odelta = 2
# 6Li
densfolder = r"/home/alexander/OneDrive/densities-6Li/1Ndensities/133MeV/60deg/"
outputfolder = (
    r"/home/alexander/Dropbox/PionPhotoProduction/results-6Li/133MeV/1bod/above/"
)

ro.main(densfolder, outputfolder)

# 3He
densfolder = r"/home/alexander/OneDrive/densities-3He/1Ndensities/132MeV/60deg/"
outputfolder = r"/home/alexander/Dropbox/PionPhotoProduction/results-3He/1bod/above/"

ro.main(densfolder, outputfolder)
# 3H
densfolder = r"/home/alexander/OneDrive/densities-3H/1Ndensities/132MeV/60deg/"

outputfolder = (
    r"/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/1bod/above/"
)
ro.main(densfolder, outputfolder)

# 4He
densfolder = r"/home/alexander/OneDrive/densities-4He/1Ndensities/133MeV/60deg/"
outputfolder = (
    r"/home/alexander/Dropbox/PionPhotoProduction/results-4He/133MeV/1bod/above/"
)

ro.main(densfolder, outputfolder)
