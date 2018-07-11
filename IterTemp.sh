#! /bin/bash

echo "start Iter Temp Fit"


rm -r IterationProgress
rm -r MCTemplatesAnData
rm -r MixedBGComp
mkdir IterationProgress
mkdir MCTemplatesAnData
mkdir MixedBGComp


time root -l -b -q IterTempPlot.C++\($1,\"$2\"\)
