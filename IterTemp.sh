#! /bin/bash

echo "start Iter Temp Fit"

# rm -r MCTemplatesAnData
# rm -r MixedBGComp
# mkdir MCTemplatesAnData
# mkdir MixedBGComp

DIR=${PWD##*/}

time root -l -b -q IterTempCreation.C++\(\"$DIR\"\)
time root -l -b -q IterTempPlot.C++\($1,\"$2\",\"$3\"\)
