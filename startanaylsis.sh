#! /bin/bash

echo ""
echo "start Template_CAP.C"

DIR=${PWD##*/} # the directery where I am currently to make things a bit more
               # felxible in the code
#$1            Which method should be used:
#$1            == 1 uses the backgrund templates from the 3 to 8 method
#              == 2 uses the  -||-                        Next Neighbours method
#
#

time root -l -q -b Template_CAP.C++\(\"$DIR\",$0\)
