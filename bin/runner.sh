#!/bin/bash
#PBS -N graphgen
#PBS -l select=1:ncpus=24:mem=500gb,walltime=72:00:00
#PBS -q bigmem
#PBS -o /home/jsybran/job.out
#PBS -e /home/jsybran/job.err
#PBS -M jsybran@clemson.edu
#PBS -m ea

LOG_FILE=/scratch2/jsybran/log.txt

export LD_LIBRARY_PATH=~/lib

module load gcc openmpi/1.10.3

rm $LOG_FILE
touch $LOG_FILE

$tail -f $LOG_FILE &
/home/jsybran/projects/MedlineGraphGenerator/bin/medlineGraphGen

