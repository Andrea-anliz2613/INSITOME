#!/bin/bash

echo "Converting ANCESTRYMAP to EIGENSTRAT"
../bin/convertf -p ../data/parfiles/par.ANCESTRYMAP.EIGENSTRAT

echo "Converting EIGENSTRAT to PED"
../bin/convertf -p ../data/parfiles/par.EIGENSTRAT.PED

echo "Finished converting to PED"
