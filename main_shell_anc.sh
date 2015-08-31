#!/bin/sh

#PBS -N MAIN_SHELL
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=1
#PBS -q short
#PBS -o outputmain_anc.file
#PBS -e errormain_anc.file

speed=(0.5 1 1.5 2)
variance=(0 2 4 6 8 10)


for S in "${speed[@]}"
do
for V in "${variance[@]}"
do
export SPE=$S
export VAR=$V
qsub sub_shell_anc.sh -v SPE,VAR
done
done
