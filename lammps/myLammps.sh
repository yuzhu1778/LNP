#!/bin/bash


python3 calCGMass.py 14734 LNPmass.dat

python3 LammpsG.py etest.xyz LNPmass.dat myLNP.dat myLNP.in
