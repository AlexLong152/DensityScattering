# no arguments ---------------------------- runs input.dat
# argument "test", or "t", just passes ---- runs testInput.dat
# argument "physical" or "p" -------------- runs physicalInput.dat
#clear
make clean
make
chmod +x run.twobodyvia2Ndensity.PionPhotoProdThresh

# if no arguments then run input.dat
# ./run.twobodyvia2Ndensity.PionPhotoProdThresh input.dat
# ./run.twobodyvia2Ndensity.PionPhotoProdThresh tIn.dat
#
if [ -z "$1" ]; then
  # No argument provided, use default input.dat
  ./run.twobodyvia2Ndensity.PionPhotoProdThresh input.dat
else
  # Argument provided, use input-<arg>.dat
  ./run.twobodyvia2Ndensity.PionPhotoProdThresh input-$1.dat
fi
