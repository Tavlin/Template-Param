#! /bin/bash

echo ""
echo "start Template_CAP.C"

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code


for i in {1..3..1}
  do

    echo ""
    echo "Start templatemethod number $i analysis"
    echo ""
    time root -l -q -b Template_CAP.C++\(\"$DIR\",$i\)
  done

# rm -r BetterBkgNN
# mkdir BetterBkgNN
# rm -r BetterBkg3to8
# mkdir BetterBkg3to8
# rm -r BetterBkg3to8Pulse
# mkdir BetterBkg3to8Pulse

# time root -l -q -b TemplatePlotting.C++\(\"$1\",\"$2\"\)
