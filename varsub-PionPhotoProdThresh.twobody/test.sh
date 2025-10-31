# no arguments ---------------------------- runs input.dat
# argument "test", or "t", just passes ---- runs testInput.dat
# argument "physical" or "p" -------------- runs physicalInput.dat
#clear
make clean
make
chmod +x run.twobodyvia2Ndensity.PionPhotoProdThresh

./run.twobodyvia2Ndensity.PionPhotoProdThresh tIn.dat
