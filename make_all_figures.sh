# make_all_figures.sh
# not actually intended to be run as-is; just to document the/an appropriate
# order in which to call all the scripts necessary to create the figures and
# which arguments I used

# only really need 4 threads to test iML1515 and yeast-GEM in a decent amount
# of time, but definitely need 8 or more for Human-GEM, and as many as you can
# get to test all the Mendoza et al. 2019 models or the AGORA2 models
threads=4

# set up the virtual environment
python -m virtualenv .
source bin/activate
pip install git+https://github.com/Devlin-Moyer/macaw.git@main

# download all of the AGORA2 models and extract the zip archives they come in
mkdir -p GSMMs/AGORA2
cd GSMMs/AGORA2
wget https://www.vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/zipped/AGORA2_annotatedMat_A_C.zip.zip
wget https://www.vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/zipped/AGORA2_annotatedMat_D_F.zip.zip
wget https://www.vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/zipped/AGORA2_annotatedMat_G_P.zip.zip
wget https://www.vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/zipped/AGORA2_annotatedMat_R_Y.zip
unzip AGORA2_annotatedMat_A_C.zip.zip
unzip AGORA2_annotatedMat_G_P.zip.zip
unzip AGORA2_annotatedMat_D_F.zip.zip
unzip AGORA2_annotatedMat_R_Y.zip
rm AGORA2_annotatedMAT_*.zip
cd ../..

# create all the data needed to make figures
mkdir -p logs
python scripts/run_tests_on_human.py 1.15 $threads > logs/Human-GEMv1.15_log.txt
python scripts/run_tests_on_human.py 1.19 $threads > logs/Human-GEMv1.19_log.txt
python scripts/run_tests_on_yeast.py 9.0.0 $threads > logs/yeast-GEMv9.0.0_log.txt
python scripts/run_tests_on_ecoli.py $threads > logs/iML1515_log.txt
python scripts/fig_2_S4b_data.py
python scripts/fig_3a_S2a_S3a_data.py
python scripts/fig_4_data.py GSMMs/mendoza_2019 fig_4a_data $threads 1 1 > logs/fig_4a_log.txt
python scripts/fig_4_data.py GSMMs/AGORA2 fig_4b_data $threads 1 1 >> logs/fig_4b_log.txt
python scripts/fig_6c_data.py $threads
python scripts/fig_S4a_data.py $threads >> logs/fig_S4a_log.txt
python scripts/fig_S5_data.py

# fig_4_data.py and fig_S4a_data.py need lots of threads (>=16) and considerable
# memory (>500GB) to finish in anything under a week

# make figures (note that all Escher maps were manually created)
# this was all done with R version 4.3.2
Rscript scripts/fig_2_S4.R
Rscript scripts/fig_3_S2_S3.R
Rscript scripts/fig_4.R
Rscript scripts/fig_6.R
Rscript -e "rmarkdown::render('scripts/additional_file_1.Rmd', output_file = '../figures/Additional File 1.pdf')"
