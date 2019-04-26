#! /bin/bash
# $a == what should be plotted
# $b == SaveFile Format (.eps, .pdf, .png etc.)
# $c == path to the ESD MC
# $d == path to the ESD data
# $e == path to the normal framework output
# $f == path to the framework output without rebinning
# $g == path to the framework output with higher rebinning
# $h == path to the framework output with lower rebinning
# $i == path to the framework output with narrow uncorr back fit range
# $j == path to the framework output with wide uncorr back fit range
# $k == cut string (since the cut string log can contain more then one and this is made only for one cut string!)
# $l == SaveAppendix (if you want to add something like the data set to the pictures' names)

source Files.txt

echo ""
echo "start Template_CAP.C"

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code


## changed loop from 1-5 to 2-5 since! needs to be changed back!

for iTM in {0..0..-1}
  do

    echo ""
    echo "start templatemethod number $iTM analysis"
    echo ""
    time root -l -q -b Template_CAP.C++\($iTM,\"$c\",\"$d\",\"$e\",\"$f\",\"$g\",\"$h\",\"$i\",\"$j\",\"$k\"\)
  done

# rm -r BetterBkgNN
# mkdir BetterBkgNN
# rm -r BetterBkg3to8
# mkdir BetterBkg3to8
# rm -r Normal
# mkdir Normal
#
# time root -l -q -b TemplatePlotting.C++\(\"$a\",\"$b\"\,\"$e\",\"$l\"\)

rm -r Yields
mkdir Yields
time root -l -q -b Yields.C++\(\"$b\",\"$l\"\)
time root -l -q -b YieldPlotting.C++\(\"$a\",\"$b\"\,\"$e\",\"$l\"\)


# rm -r Untergrund
# mkdir Untergrund
# time root -l -q -b CorrBackFitPlot.C++\(\"$b\"\)
