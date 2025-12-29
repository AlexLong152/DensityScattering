# Compile gwwi08.f with legacy Fortran flags
rm piScat
gfortran -std=legacy -fallow-argument-mismatch -o piScat gwwi08.f
./piScat
