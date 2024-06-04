#!/bin/bash -l
#$ -m e
#$ -j y
#$ -l h_rt=24:00:00
#$ -pe omp 16
#$ -t 1-2

source bin/activate

if [ $SGE_TASK_ID == 1 ] 
then
  python run_tests_on_human.py 1.15 8 > logs/Human-GEMv1.15_log.txt
else
  python run_tests_on_human.py 1.18 8 > logs/Human-GEMv1.18_log.txt
fi
