#! /bin/bash

echo ""
echo "start Template_CAP.C"

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code


for i in {1..4..1}
  do

    echo ""
    echo "Start templatemethod number $i analysis"
    echo ""
    time root -l -q -b Template_CAP.C++\(\"$DIR\",$i\)
  done

time root -l -q -b Template_CAP.C++\(\"$DIR\",4\)

rm -r BetterBkgNN
mkdir BetterBkgNN
rm -r BetterBkg3to8
mkdir BetterBkg3to8
rm -r BetterBkg3to8Pulse
mkdir BetterBkg3to8Pulse
rm -r Normal
mkdir Normal

time root -l -q -b TemplatePlotting.C++\(\"$1\",\"$2\"\)

rm -r Systematics
mkdir Systematics
time root -l -q -b Systematics.C++\(\"$2\"\)

# rm -r Untergrund
# mkdir Untergrund
# time root -l -q -b CorrBackFitPlot.C++\(\"$2\"\)
