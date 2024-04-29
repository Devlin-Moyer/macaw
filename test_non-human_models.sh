#!/bin/bash -l
#$ -m e
#$ -j y
#$ -l h_rt=12:00:00
#$ -pe omp 16

source bin/activate

python run_tests_on_ecoli.py $NSLOTS > figure_data/iML1515_log.txt
python run_tests_on_yeast.py 9.0.0 $NSLOTS > figure_data/yeast-GEMv9.0.0_log.txt
