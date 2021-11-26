#!/bin/bash
#$ -cwd
#$ -V -S /bin/bash
#$ -q mag2
#$ -N testjob
#$ -pe smp 48
#$ -o kkr.log
#$ -e kkr.err

ulimit -s unlimited
#export KMP_STACKSIZE=2g

export OMP_NUM_THREADS=$NSLOTS

/usr/local/anaconda3/bin/python -u ./testrun1.py > std1.out 2>&1
/usr/local/anaconda3/bin/python -u ./testrun2.py > std2.out 2>&1
