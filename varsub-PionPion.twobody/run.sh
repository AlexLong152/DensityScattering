# no arguments ---------------------------- runs input.dat
# argument "test", or "t", just passes ---- runs testInput.dat
# argument "physical" or "p" -------------- runs physicalInput.dat
#clear
make clean
make
chmod +x run.twobodyvia2Ndensity.PionPion

# if no arguments then run input.dat
./run.twobodyvia2Ndensity.PionPion input.dat
