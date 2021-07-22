#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere 
# else in the code except when reading the input options.
#############################################################################

import sys
import timeit
import lib.core as PC

def init():
    globs = {
        'phyloacc-versions' : '',
        'interface-version' : 'Beta 1.0',
        'releasedate' : 'NA',
        'doi' : 'https://doi.org/10.1093/molbev/msz049',
        'http' : 'https://xyz111131.github.io/PhyloAcc/',
        'github' : 'https://github.com/xyz111131/PhyloAcc/issues',
        'starttime' : timeit.default_timer(),
        'startdatetime' : PC.getOutTime(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        # Script call info

        'aln-file' : False,
        'bed-file' : False,
        'id-file' : False,
        # Input with concatenated alignment and partitions by bed file

        'aln-dir' : False,
        # Input by a directory full of individual alignments

        'mod-file' : False,
        # Input rate and tree file from PHAST

        'outdir' : '',
        'run-name' : 'phyloacc',
        'logfilename' : 'phyloacc.errlog',
        'logdir' : '',
        'tmpdir' : 'System default.',
        'overwrite' : False,
        # I/O options
        
        'targets' : False,
        'conserved' : False,
        'outgroup' : False,
        # Phylo options

        'burnin' : 500,
        'mcmc' : 1000,
        'chain' : 1,
        # MCMC options

        'phyloacc' : 'PhyloAcc/PhyloAcc',
        'phyloacc-gbgc' : 'PhyloAcc/V2_GBGC/PhyloAcc_gBGC',
        'phyloacc-gt' : '',
        # Dependency paths

        'num-jobs' : 1,
        # Number of jobs/threads to use

        'norun' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
        'log-v' : 1,
        'stats' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : "",
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs['logfilename'] = "phyloacc-" + globs['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    return globs;