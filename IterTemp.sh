#! /bin/bash

echo "start Iter Temp Fit"

# rm -r MCTemplatesAnData
# rm -r MixedBGComp
mkdir MCTemplatesAnData
mkdir BackGroundFitting
mkdir Systematics
# mkdir MixedBGComp
# rm -r DPG
mkdir DPG
mkdir Variation
mkdir $5

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code

# $1 give existing bin to only plot plots for this bin, otherwise, all bins.
# $2 string to call which of the plots should be done. all for all obv.
# $3 save format. png or eps
# $4 mode 1 for p-Pb else for pp
# $5 Saving in another Diretory, for like testing and stuff.

# time root -l -b -q IterTempCreation.C++\(\"$DIR\",$4\)
# time root -l -b -q IterTempPlot.C++\($1,\"$2\",\"$3\",\"$5\"\)
# time root -l -b -q dpg.C++\($1,\"$2\",\"$3\",\"$5\"\)
# time root -l -b -q BackGroundFitting.C++\(\"$3\"\)
time root -l -b -q Systematics.C++\(\"$3\"\)
