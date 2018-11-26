#! /bin/bash
# $a == what should be plotted
# $b == SaveFile Format (.eps, .pdf, .png etc.)
# $c == ESD_MC
# $d == ESD_data
# $e == MC
# $f == Data
# $g == CorrectedData
# $h == Correction
# $x == MCRebin1
# $SaveAppendix

source Files.txt

echo ""
echo "start Template_CAP.C"

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code


## changed loop from 1-5 to 2-5 since! needs to be changed back!

for i in {1..5..1}
  do

    echo ""
    echo "Start templatemethod number $i analysis"
    echo ""
    time root -l -q -b Template_CAP.C++\(\"$DIR\",$i,\"$c\",\"$d\",\"$e\",\"$f\",\"$g\",\"$h\",\"$x\"\)
  done

# time root -l -q -b Template_CAP.C++\(\"$DIR\",3,\"$c\",\"$d\",\"$e\",\"$f\,\"$g\",\"$h\",\"$x\"\)

# rm -r BetterBkgNN
# mkdir BetterBkgNN
# rm -r BetterBkg3to8
# mkdir BetterBkg3to8
# rm -r BetterBkg3to8Pulse
# mkdir BetterBkg3to8Pulse
# rm -r Normal
# mkdir Normal
# rm -r OneTemplate
# mkdir OneTemplate
#
# time root -l -q -b TemplatePlotting.C++\(\"$a\",\"$b\"\,\"$f\",\"$SaveAppendix\"\)

rm -r Yields
mkdir Yields
time root -l -q -b Yields.C++\(\"$b\",\"$SaveAppendix\"\)


# rm -r Untergrund
# mkdir Untergrund
# time root -l -q -b CorrBackFitPlot.C++\(\"$b\"\)
