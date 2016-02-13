#!/bin/bash 

for bctyp in {bc,bcm,bch,bchh}
do
    echo $bctyp
    python2 cmpr_mec.py $bctyp
done
