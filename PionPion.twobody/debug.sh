rm *.tmp
clear;
make clean
cp Makefile makedebug
#replace -O optimizeation with no optimization at all for better debugging
#replaces the string "FFLAGS= -O -" with the string "FFLAGS= -O0 -"
sed -i 's/FFLAGS= -O -/FFLAGS= -O0 -/' makedebug
make -f makedebug # runs make with the debug make file
# now remove the makedebug file to avoid confusion
rm makedebug
chmod +x run.twobodyvia2Ndensity
echo ""
echo "Debuging with gdb debugger"
echo "##################################################################"
echo "Set breakpoints with------------------b filename.f:linenumber"
echo "Run code until breakpoint with--------r or run"
echo "Run until next breakpoint with--------c or continue"
echo "Print with----------------------------p:variablename"
echo "Go to next line in file with----------n"
echo "Step into function with---------------s or step"
echo ""
echo "More infor here: https://cs.baylor.edu/~donahoo/tools/gdb/tutorial.html"
echo ""
gdb --args twobodyvia2Ndensity testInput.dat
