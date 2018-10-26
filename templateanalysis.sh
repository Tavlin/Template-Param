#! /bin/bash

echo ""
echo "start Template_CAP.C"

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code

time root -l -q -b Template_CAP.C++\(\"$DIR\",$2\)
