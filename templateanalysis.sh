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

# time root -l -q -b Template_CAP.C++\(\"$DIR\",3\)
