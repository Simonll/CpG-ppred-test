#!/bin/bash

echo $1 

rep=$1-$(date +"%D"-"%T"| sed "s/\//-/g"| sed "s/:/-/g")
cd /tmp/
git clone --recursive git@github.com:Simonll/$1.git $rep
zip -r $rep.zip $rep -x '*.git*'  $rep/src/\* $rep/notes/\* $rep/notebooks/\* $rep/data/\*
