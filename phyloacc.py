#!/usr/bin/env python3
#############################################################################
# This is the python front-end for PhyloAcc, a Bayesian substitution rate
# estimation program for conserved non-coding genomic elements. This script
# will handle user inputs, model selection, and batching jobs
#
# Gregg Thomas
# Summer 2021
#############################################################################

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

#############################################################################

import sys
import os
import lib.core as PC
import lib.params as params
import lib.opt_parse as OP
import lib.seq as SEQ
import lib.tree as TREE
import lib.output as OUT
import lib.batch as BATCH

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.
    
    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# PhyloAcc version " + globs['interface-version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual PhyloAcc version for this, and not just the interface version.

    print("#");
    print("# " + "=" * 125);
    print(PC.welcome());
    if "-h" not in sys.argv:
        print("    Bayesian rate analysis of conserved");
        print("       non-coding genomic elements\n")
    # A welcome banner.

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    step_start_time = PC.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    globs = SEQ.readSeq(globs);
    # Library to read input sequences

    globs = SEQ.alnStats(globs);
    # Calculate some basic alignment stats

    globs = TREE.scf(globs);
    # Calculate avg. sCF per locus

    globs = OUT.writeAlnStats(globs);
    # Write out the alignment summary stats

    globs = OUT.writeSCFStats(globs);
    # Write out the sCF summary stats

    globs = BATCH.genJobFiles(globs);
    # Generates the locus specific job files (aln, bed, config, etc.) for phyloacc

    globs = BATCH.writeSnakemake(globs);
    # Generates the snakemake config and cluster profile

    PC.endProg(globs);

#############################################################################

