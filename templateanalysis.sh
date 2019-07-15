#! /bin/bash

: ' $a == what should be plotted
    $b == SaveFile Format (.eps, .pdf, .png etc.)
    $c == path to the ESD MC
    $d == path to the ESD data
    $e == path to the normal framework output
    $f == path to the framework output without rebinning
    $g == path to the framework output with higher rebinning
    $h == path to the framework output with lower rebinning
    $i == path to the framework output with narrow uncorr back fit range
    $j == path to the framework output with wide uncorr back fit range
    $k == cut string (since the cut string log can contain more then one and this is made only for one cut string!)
    $l == SaveAppendix (if you want to add something like the data set to the pictures names)
'

source Files.txt

echo ""
echo "START TEMPLATE ANALYSIS..."

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code


: ' This loops over all needed "analysis tasks". The numbers correspond to one
    task each:
              [==  0 uses the background templates from the 3 to 8 method,
               ==  1 uses the background templates from the Next Neighbours method
               ==  2 uses the background templates from the same pT intervals as the signal templates
               ==  3 uses lower rebinning output              with TEMPLATEMETHOD 1
               ==  4 uses higher rebinning output             with TEMPLATEMETHOD 1
               ==  5 uses narrow uncorr. bck fit range output with TEMPLATEMETHOD 1
               ==  6 uses wide uncorr. bck   fit range output with TEMPLATEMETHOD 1
               ==  7 uses narrow             fit range output with TEMPLATEMETHOD 1
               ==  8 uses wide               fit range output with TEMPLATEMETHOD 1
               ==  9 uses narrow           count range output with TEMPLATEMETHOD 1
               == 10 uses wide             count range output with TEMPLATEMETHOD 1]
'

for iTM in {10..0..-1}
  do

    echo "_____________________________________________________________________"
    echo "start templatemethod number $iTM analysis task..."
    echo ""
    time root -l -q -b Template_CAP.C++\($iTM,\"$c\",\"$d\",\"$e\",\"$f\",\"$g\",\"$h\",\"$i\",\"$j\",\"$k\"\)
  done

echo "_____________________________________________________________________"
echo "removing and recreating folders for monitoring plots"
echo ""

rm -r BetterBkgNN
mkdir BetterBkgNN
rm -r BetterBkg3to8
mkdir BetterBkg3to8
rm -r Normal
mkdir Normal

echo "_____________________________________________________________________"
echo "start plotting macro for general monitoring..."

time root -l -q -b TemplatePlotting.C++\(\"$a\",\"$b\",\"$e\",\"$k\",\"$l\"\)

echo "_____________________________________________________________________"
echo "removing and recreating folders for yield and systematic uncertainties
plots"
echo ""

rm -r Yields
mkdir Yields

echo "_____________________________________________________________________"
echo "start plotting macro for yield and systematic uncertainties..."

time root -l -q -b Yields.C++\(\"$b\",\"$l\"\)
time root -l -q -b YieldPlotting.C++\(\"$a\",\"$b\"\,\"$e\",\"$l\"\)
