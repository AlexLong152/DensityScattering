# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import runfolderOnebody as ro
from os import system

system("make clean;make")
ro.Odelta = 2
# 6Li
"""
densfolder = r"/home/alexander/OneDrive/densities-6Li/1Ndensities/138MeV/"

outputfolder = r"/home/alexander/Dropbox/PionPion/results-6Li/1bod/thresh"

ro.main(densfolder, outputfolder)

# 3He
densfolder = r"/home/alexander/OneDrive/densities-3He/1Ndensities/136MeV/"
outputfolder = "/home/alexander/Dropbox/PionPion/results-3He/1bod/thresh/"

ro.main(densfolder, outputfolder)

# 3H
# densfolder = r"/home/alexander/OneDrive/densities-3H/1Ndensities/132MeV/60deg/"
# outputfolder = (
#     r"/home/alexander/Dropbox/PionPhotoProduction/results-3H/132MeV/1bod/above/"
# )
# ro.main(densfolder, outputfolder)

# 4He
densfolder = r"/home/alexander/OneDrive/densities-4He/1Ndensities/136MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-4He/1bod/thresh/"
ro.main(densfolder, outputfolder)

#############################################################################################
##################      pi0 stuff now     ###################################################
#############################################################################################
# 6Li
densfolder = r"/home/alexander/OneDrive/densities-6Li/1Ndensities/133MeV/0deg/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-6Li/1bod/threshpi0/"

ro.main(densfolder, outputfolder)

# 3He
densfolder = r"/home/alexander/OneDrive/densities-3He/1Ndensities/132MeV/0deg"
outputfolder = "/home/alexander/Dropbox/PionPion/results-3He/1bod/threshpi0/"

ro.main(densfolder, outputfolder)


# 4He
densfolder = r"/home/alexander/OneDrive/densities-4He/1Ndensities/133MeV/0deg/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-4He/1bod/threshpi0/"
ro.main(densfolder, outputfolder)
"""

# 4He, above threshold
densfolder = r"/home/alexander/OneDrive/densities-4He/1Ndensities/PionScat-168MeV/"
outputfolder = r"/home/alexander/Dropbox/PionPion/results-4He/1bod/170MeV-pi0/"
ro.main(densfolder, outputfolder)
