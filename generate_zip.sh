#!/bin/bash

echo $1

rep=$1-$(date +"%D"-"%T"| sed "s/\//-/g"| sed "s/:/-/g")
cd /tmp/
git clone --recursive git@github.com:Simonll/$1.git $rep
zip -r $rep.zip $rep/main.tex $rep/figures/ $rep/sections $rep/refs.bib $rep/tables $rep/spbasic.bst $rep/main.pdf
