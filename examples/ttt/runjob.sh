#!/bin/bash
#PBS -l nodes=1:ppn=16
##PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
###PBS -q roux
#PBS -j oe
#PBS -A Drude
#PBS -N testGAAMP

cd $PBS_O_WORKDIR

export TMPDIR="/home/eliot/tmp"

cat $PBS_NODEFILE > nodes.txt

/home/eliotblg/work/GAAMP-17/scripts/gaamp ethanol.inp

