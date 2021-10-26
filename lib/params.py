#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere 
# else in the code except when reading the input options.
#############################################################################

import sys
import timeit
import lib.core as PC

#############################################################################

class StrictDict(dict):
# This prevents additional keys from being added to the global params dict in
# any other part of the code, just to help me limit it's scope
# https://stackoverflow.com/questions/32258706/how-to-prevent-key-creation-through-dkey-val
    def __setitem__(self, key, value):
        if key not in self:
            raise KeyError("{} is not a legal key of this StrictDict".format(repr(key)));
        dict.__setitem__(self, key, value);

#############################################################################

def init():
    globs_init = {
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

        'seq-compression' : 'none',
        'bed-compression' : 'none',
        # The type of compression used for input sequence files

        'outdir' : '',
        'run-name' : 'phyloacc',
        'logfilename' : 'phyloacc.errlog',
        'alnstatsfile' : 'phyloacc-aln-stats.csv',
        'scfstatsfile' : 'phyloacc-scf-stats.csv',
        'scftreefile' : 'phyloacc-scf.tree',
        'logdir' : '',
        'tmpdir' : 'System default.',
        'overwrite' : False,
        # I/O options
        
        'tree-string' : False,
        'tree-dict' : False,
        'labeled-tree' : False,
        'root-node' : False,
        'tree-tips' : False,
        # Tree variables

        'in-seqs' : {},
        'in-bed' : {},
        'alns' : {},
        'aln-stats' : {},
        'num-loci' : False,
        # Sequence variables

        'targets' : [],
        'conserved' : [],
        'outgroup' : [],
        # Phylo options

        'scf' : {},
        # Tracks sCF values BY NODE over all loci

        'quartets' : {},
        # All qurtets for all nodes/branches in the species tree:
        # <node> : [ <quartets> ]
        # <quartet> is a tuple of tuples:
        # ((split1-spec1, split1-spec2), (split2-spec1, split2-spec2))

        'min-scf' : 0.5,
        # sCF params

        'burnin' : 500,
        'mcmc' : 1000,
        'chain' : 1,
        # MCMC options

        'phyloacc' : 'PhyloAcc/PhyloAcc',
        'phyloacc-gbgc' : 'PhyloAcc/V2_GBGC/PhyloAcc_gBGC',
        'phyloacc-gt' : 'PhyloAcc-GT2/SRC/PhyloAcc-GT_piQ',
        # Dependency paths

        'batch-size' : 50,
        'num-batches' : 0,
        # Batch variables

        'num-procs' : 1,
        # Number of procs for this script to use

        'phyloacc-procs' : 1,
        'num-jobs' : 1,
        'procs-per-job' : 1,
        # Number of jobs/procs for PhyloAcc to use

        'partition' : False,
        'num-nodes' : "1",
        'mem' : "4",
        'time' : "1:00:00",
        # Cluster options

        'aln-pool' : False,
        'scf-pool' : False,
        # Process pools

        'smk' : False,
        'smk-config' : False,
        # Job files

        'job-dir' : '',
        'job-alns' : '',
        'job-cfgs' : '',
        'job-bed' : '',
        'job-smk' : '',
        'job-out' : '',
        'profile-dir' : False,
        # Job directories

        'label-tree' : False,
        'info' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

        'skip-chars' : ["-", "N"],
        'aln-stats-written' : False,
        'scf-stats-written' : False,
        'scf-tree-written' : False,
        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
        'log-v' : 1,
        'stats' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : False,
        'qstats' : False,
        'norun' : False,
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs_init['logfilename'] = "phyloacc-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################