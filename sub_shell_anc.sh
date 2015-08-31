#!/bin/sh

#PBS -N SUB_SHELL
#PBS -l walltime=02:59:00
#PBS -l nodes=1:ppn=1
#PBS -q short
#PBS -e suberror_anc.file

module load scripts
module load Python/2.7.2-ictce-4.0.6
chmod 770 anc_model.py
python anc_model.py $SPE $VAR
