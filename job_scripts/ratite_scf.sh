#!/bin/bash
#SBATCH --job-name=phyloacc_ratite_scf
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gthomas@g.harvard.edu
#SBATCH --partition=holy-info,holy-smokes,holy-cow
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=12g
#SBATCH --time=12:00:00

padir="/n/home07/gthomas/env/pkgs/PhyloAcc-interface/"
alnfile="data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz"
bedfile="data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed.gz"
modfile="data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod"
idfile="data/ratite_data/07_cnees/datasets/original_dataset_v2/ratited-ids-2000.txt"
targets="strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid"
outgroups="allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar"

cd $padir

phyloacc.py -a $alnfile -b $bedfile -m $modfile -i $idfile -o test-real-2000 -t "$targets" -g "$outgroups" -p 32 --overwrite