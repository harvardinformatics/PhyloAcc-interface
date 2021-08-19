#!/usr/bin/env python3
#############################################################################
# This is the python front-end for PhyloAcc, a Bayesian substitution rate
# estimation program for conserved non-coding genomic elements. This script
# will handle user inputs, model selection, and batching jobs
#
# Gregg Thomas
# Summer 2021
#############################################################################

# phyloacc.py -a PhyloAcc/Simulation_ratite/simu_500_200_diffr_2-1.fasta -b PhyloAcc/Simulation_ratite/simu_500_200_diffr_2-1.bed -m PhyloAcc/Data/ratite/neut_ver3_final.named.mod -o test -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 12 --overwrite

# phyloacc.py -a data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concatenated_v2.fasta.gz -b data/ratite_data/07_cnees/datasets/original_dataset_v2/allspecies_cnee_concat_partitions.bed.gz -m data/ratite_data/07_cnees/datasets/original_dataset_v2/neut_final_orig_v2.named.mod -o test-real -t "strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid" -g "allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar" -p 12 --overwrite
 
import sys
import os
import lib.core as PC
import lib.params as params
import lib.opt_parse as OP
import lib.seq as SEQ
import lib.tree as TREE
import lib.output as OUT

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

    PC.endProg(globs);

#############################################################################

