# no arguments ---------------------------- runs input.dat
# argument "test", or "t", just passes ---- runs testInput.dat
# argument "physical" or "p" -------------- runs physicalInput.dat
#clear
make clean
make
chmod +x run.twobodyvia2Ndensity.PionPhotoProdThresh

if [ $# -eq 0 ]; then
# if no arguments then run input.dat
./run.twobodyvia2Ndensity.PionPhotoProdThresh input.dat
fi

if [[ "$1" == "test" || "$1" == "t" ]]
then # if "test" is passed then run testInput.dat
./run.twobodyvia2Ndensity.PionPhotoProdThresh testInput.dat
echo ""
echo "Compare the above to GOSPEL trivial integration"
sed -n '726,737p;738q' results/GOSPEL-mpi-trivial.dat
fi

if [[ "$1" == "p" || "$1" == "physical" ]]
then # if "physical" is passed then run physicalInput.dat
./run.twobodyvia2Ndensity.PionPhotoProdThresh physicalInput.dat
fi
