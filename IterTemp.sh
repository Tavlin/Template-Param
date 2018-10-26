#! /bin/bash

echo "start Iter Temp Fit"

# rm -r MCTemplatesAnData
# mkdir MCTemplatesAnData
# rm -r MixedBGComp
# mkdir MixedBGComp
# rm -r BackGroundFitting
# mkdir BackGroundFitting
# rm -r Systematics
# mkdir Systematics
# rm -r DPG
# mkdir DPG
# rm -r Variation
# mkdir Variation
rm -r $5
mkdir $5

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code

# $1 give existing bin to only plot plots for this bin, otherwise, all bins.
# $2 string to call which of the plots should be done. all for all obv.
# $3 save format. png or eps
# $4 mode 1 for p-Pb else for pp
# $5 Saving in another Diretory, for like testing and stuff.

# for i in {2..38..2}
#   do
#   time root -l -b -q IterTempCreation2.C++\(\"$DIR\",$4,$i\)
#   # time root -l -b -q IterTempPlot.C++\($1,\"$2\",\"$3\",\"$5\"\)
#   # time root -l -b -q dpg.C++\($1,\"$2\",\"$3\",\"$5\"\)
#   # time root -l -b -q BackGroundFitting.C++\(\"$3\"\)
#   # time root -l -b -q BkgPlotting.C++\(\"$3\"\)
#   time root -l -b -q Systematics.C++\(\"$3\",$i\)
#  done

time root -l -b -q IterTempCreation2.C++\(\"$DIR\",$4,6\)

# rm -r Systematics
# mkdir Systematics
# time root -l -b -q Systematics.C++\(\"$3\",6\)
#
# rm -r BetterBkgNN
# mkdir BetterBkgNN
# rm -r BetterBkg3to8
# mkdir BetterBkg3to8
# rm -r Normal
# mkdir Normal
# time root -l -b -q IterTempPlot.C++\($1,\"$2\",\"$3\",\"$5\"\)

# time root -l -b -q BackGroundFitting.C++\(\"$3\"\)
# time root -l -b -q BkgPlotting.C++\(\"$3\"\)

# time root -l -b -q YieldStatUnc.C++\(\"$3\"\)

# time root -l -b -q MiniSim.C++\(\)

# time root -l -b -q Ratio.C++\(\)
