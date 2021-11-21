#!/usr/bin/env python3
#############################################################################
# This is the python front-end for PhyloAcc, a Bayesian substitution rate
# estimation program for conserved non-coding genomic elements. This script
# will handle user inputs, model selection, and batching jobs
#
# Gregg Thomas
# Summer 2021
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
    
    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# PhyloAcc version " + globs['interface-version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual PhyloAcc version for this, and not just the interface version.

    if "--quiet" not in sys.argv:
        print("\n" + " ".join(sys.argv) + "\n");
        # Print the command used to call the program

        print("#");
        print("# " + "=" * 125);
        print(PC.welcome());
        if "-h" not in sys.argv:
            print("    Bayesian rate analysis of conserved");
            print("       non-coding genomic elements\n");
        # A welcome banner.
    else:
        print("# Running phyloacc_interface in quiet mode...");

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

    if globs['run-mode'] == 'adaptive':
        globs = TREE.scf(globs);
    # Calculate avg. sCF per locus

    globs = OUT.writeAlnStats(globs);
    # Write out the alignment summary stats

    if globs['run-mode'] == 'adaptive':
        globs = OUT.writeSCFStats(globs);
    # Write out the sCF summary stats

    globs = BATCH.genJobFiles(globs);
    # Generates the locus specific job files (aln, bed, config, etc.) for phyloacc

    globs = BATCH.writeSnakemake(globs);
    # Generates the snakemake config and cluster profile

    PC.endProg(globs);

#############################################################################

