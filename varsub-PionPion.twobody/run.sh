#clear
make clean
make
chmod +x run.twobodyvia2Ndensity.PionPion

if [ -z "$1" ]; then
  # No argument provided, use default input.dat
  ./run.twobodyvia2Ndensity.PionPion input.dat
else
  # Argument provided, use input-<arg>.dat
  ./run.twobodyvia2Ndensity.PionPion input-$1.dat
fi
