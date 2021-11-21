#!/bin/bash

time -p python phyloacc_interface.py -d data/simu_500_200_diffr_2-1/ -m PhyloAcc/Data/ratite/neut_ver3_final.named.mod -o test-adaptive-theta -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -part "holy-info,holy-cow,holy-smokes" -n 24 -p 4 -j 12 --overwrite --theta -r adaptive
time -p python phyloacc_interface.py -d data/simu_500_200_diffr_2-1/ -m PhyloAcc/Data/ratite/neut_ver3_final.named.mod -o test-adaptive -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -part "holy-info,holy-cow,holy-smokes" -n 24 -p 4 -j 12 -l sim-tree-astral.tre --overwrite -r adaptive
time -p python phyloacc_interface.py -d data/simu_500_200_diffr_2-1/ -m PhyloAcc/Data/ratite/neut_ver3_final.named.mod -o test-adaptive -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -part "holy-info,holy-cow,holy-smokes" -n 24 -p 4 -j 12 --overwrite
time -p python phyloacc_interface.py -d data/simu_500_200_diffr_2-1/ -m PhyloAcc/Data/ratite/neut_ver3_final.named.mod -o test-gt -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -part "holy-info,holy-cow,holy-smokes" -n 24 -p 4 -j 12 -l sim-tree-astral.tre --overwrite -r gt
time -p python phyloacc_interface.py -d data/simu_500_200_diffr_2-1/ -m PhyloAcc/Data/ratite/neut_ver3_final.named.mod -o test-gt-theta -t strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid -g allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar -part holy-info,holy-cow,holy-smokes -n 24 -p 4 -j 12 --theta --overwrite -r gt

##########################
## OLD CMDS
##########################

# Simulated dataset:

# time -p python phyloacc.py -d data/simu_500_200_diffr_2-1/ -m PhyloAcc/Data/ratite/neut_ver3_final.named.mod -o test -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -part "holy-info,holy-cow,holy-smokes" -n 24 -p 24 -j 6
# real 70.09
# user 1608.92
# sys 3.24

# time -p snakemake -p -s test/phyloacc-job-files/snakemake/run_phyloacc.smk --configfile test/phyloacc-job-files/snakemake/phyloacc-config.yaml --profile test/phyloacc-job-files/snakemake/profiles/slurm_profile --dryrun
# real 4917.29
# user 1.82
# sys 0.60

##########################

# Simulated dataset (MSC):

# time -p python phyloacc.py -a PhyloAcc-GT2/simu_MSC/simu_200_100_2-5.fasta -b PhyloAcc-GT2/simu_MSC/simu_200_100_2-5.bed -m PhyloAcc-GT2/simu_MSC/neut_ver3_genetree2.named.mod -o test-msc -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "galGal;nipNip;halLeu" -part "holy-info,holy-cow,holy-smokes" -n 24 -p 24 -j 6
# real 3.64
# user 63.67
# sys 0.61

# time -p snakemake -p -s test-msc/phyloacc-job-files/snakemake/run_phyloacc.smk --configfile test-msc/phyloacc-job-files/snakemake/phyloacc-config.yaml --profile test-msc/phyloacc-job-files/snakemake/profiles/slurm_profil

##########################

# Real dataset:

# 1000 loci:
# time -p python phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -i data/ratite_data/07_cnees/datasets/original_dataset_v2/ratited-ids-1000.txt -o test-real-1000 -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 32 -part 'holy-info,holy-cow,holy-smokes,shared' --overwrite
# real 206.86
# user 2335.62
# sys 15.92

# 2000 loci
# time -p python phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -i data/ratite_data/07_cnees/datasets/original_dataset_v2/ratited-ids-2000.txt -o test-real-2000 -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 32 -part 'holy-info,holy-cow,holy-smokes,shared' --overwrite
# real 290.53
# user 4277.12
# sys 25.56

# 5000 loci
# time -p python phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -i data/ratite_data/07_cnees/datasets/original_dataset_v2/ratite-ids-5000.txt -o test-real-5000 -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 32 -part 'holy-info,holy-cow,holy-smokes,shared' --overwrite
# real 577.40
# user 10476.72
# sys 54.06

# 10000 loci
# time -p python phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -i data/ratite_data/07_cnees/datasets/original_dataset_v2/ratited-ids-10000.txt -o test-real-10000 -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 32 -part 'holy-info,holy-cow,holy-smokes,shared' --overwrite
# real 1072.99
# user 20932.03
# sys 103.50

# 50000 loci
# time -p python phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -i data/ratite_data/07_cnees/datasets/original_dataset_v2/ratite-ids-50000.txt -o test-real-50000 -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 32 -part 'holy-info,holy-cow,holy-smokes,shared' --overwrite
# real 5770.43
# user 117375.58
# sys 921.94

# 100000 loci
# time -p python phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -i data/ratite_data/07_cnees/datasets/original_dataset_v2/ratite-ids-100000.txt -o test-real-100000 -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 32 -part 'holy-info,holy-cow,holy-smokes,shared' --overwrite
# real 10343.17
# user 228212.60
# sys 1788.44

# ALL loci
# python phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed.gz -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -o /n/holylfs05/LABS/informatics/Users/gthomas/PhyloAcc-interface-data/test-real-all-b10000/ -t strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid -g allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar -n 48 -p 48 -j 48 -batch 10000 -part holy-info,holy-cow,holy-smokes,shared -time 4 --overwrite
# Total execution time:            22006.353 seconds.