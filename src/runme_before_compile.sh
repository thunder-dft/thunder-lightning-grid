#!/bin/bash

MASTER=fireball-master

for i in Makefile MACHINES j.ASSEMBLERS e.FDATA g.XC_FUNCTIONALS h.SOLVESH l.SCF o.OUTPUT u.UTIL libs include p.HARRIS d.FUNCTIONS_EXTRA a.GLOBAL b.FUNCTIONS c.SYSTEM
do
    if [ -e $i ]
    then
	rm $i
    fi
    ln -s ../../${MASTER}/src/$i
done
