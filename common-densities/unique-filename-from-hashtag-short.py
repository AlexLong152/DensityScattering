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
densdb=access.database(workdir=os.environ["HOME"]+"/EFFNUCLEON/COMPTON/FEW-NUCLEON/DENSITIES/",webbase="https://datapub.fz-juelich.de/anogga/files/densindx/")

# print all entries but restrict tto limited number of columns

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
