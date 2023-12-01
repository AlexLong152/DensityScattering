# hgrie Nov 2023:, based on files by Alexander Long and Andreas Nogga
# get the uniquefilename from a specific hashtag
# needs the indiex file densities_table.h5 of densities
# or an internet connection to download it into directory common-densities/
#
# This construcion should *NOT* need any fiddling with individual directory structures.
#
# run via >> python3.8 find_uniquefilename_from_hashtag.py <hashtag>
#
# Is also called inside fortran density code by read2Ndensity()
# to provide the unique filename to stdout if hashtag provided as part of density filename.
#
import sys
import os 
import pandas as pd

class DevNull:
    def write(self, msg):
        pass

sys.stderr = DevNull()

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

from nucdens import access


# open database using the mirror on datapub  
densdb=access.database(workdir="./",webbase="https://datapub.fz-juelich.de/anogga/files/densindx/") # original by Andreas

# print all entries but restrict to limited number of columns

colsel=["kind","Z","N","MODENN","orderNN","LambdaNN","tnforder","c1","c3","c4","cd","ce","cE1","cE3","lambdaSRGNN","srgwf","relcalc","omega","theta","qmom","uniquefilename"]

# use stdin if it's full                                                        
if not sys.stdin.isatty():
    input_stream = sys.stdin
    
# otherwise, read the given filenam
else:
    try:
        hashtag = sys.argv[1]
    except IndexError:
        message = 'need filename as first argument if stdin is not full'
        raise IndexError(message)

# then finally, do the opposite: get infos for a specific hashname
# there is only one entry possible 
# directly transfer to a dictionary
properties_dict=densdb.pddf[densdb.pddf.hashname==hashtag][colsel].to_dict('records')
# this is all info 
#print(properties_dict[0])
# but one could also look at specific entries
#print("omega,theta: ", properties_dict[0]["omega"],properties_dict[0]["theta"])
print(properties_dict[0]["uniquefilename"])
