#!/bin/bash -l
#$ -m e
#$ -j y
#$ -l h_rt=12:00:00
#$ -pe omp 8

source bin/activate
module load R

python fig_3c_data.py $NSLOTS
Rscript fig_3c.R
