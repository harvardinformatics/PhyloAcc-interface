#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import math
import argparse
import multiprocessing as mp
import lib.core as PC
import lib.tree as TREE

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py

    try:
        import psutil
        globs['psutil'] = True;
    except:
        globs['psutil'] = False;
    # Check if psutil is installed for memory usage stats.

    parser = argparse.ArgumentParser(description="PhyloAcc: Bayesian rate analysis of conserved non-coding genomic elements");

    parser.add_argument("-a", dest="aln_file", help="An alignment file with all loci concatenated. -b must also be specified. Expected as FASTA format for now. One of -a/-b or -d is REQUIRED.", default=False);
    parser.add_argument("-b", dest="bed_file", help="A bed file with coordinates for the loci in the concatenated alignment file. -a must also be specified. One of -a/-b or -d is REQUIRED.", default=False);
    parser.add_argument("-i", dest="id_file", help="A text file with locus names, one per line, corresponding to regions in the input bed file. -a and -b must also be specified.", default=False);

    parser.add_argument("-d", dest="aln_dir", help="A directory containing individual alignment files for each locus. Expected as FASTA format for now. One of -a/-b or -d is REQUIRED.", default=False);
    
    parser.add_argument("-m", dest="mod_file", help="A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from PHAST. REQUIRED.", default=False);
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: phyloacc-[date]-[time]", default=False);
    # Output
    
    parser.add_argument("-t", dest="targets", help="Tip labels in the input tree to be used as target species. Enter multiple labels separated by semi-colons (;). REQUIRED.", default=False);
    parser.add_argument("-c", dest="conserved", help="Tip labels in the input tree to be used as conserved species. Enter multiple labels separated by semi-colons (;). Any species not specified in -t or -g will be inferred as conserved.", default=False);
    parser.add_argument("-g", dest="outgroup", help="Tip labels in the input tree to be used as outgroup species. Enter multiple labels separated by semi-colons (;).", default=False);
    # Phylo options 
    
    parser.add_argument("-burnin", dest="burnin", help="The number of steps to be discarded in the Markov chain as burnin. Default: 500", default=False);
    parser.add_argument("-mcmc", dest="mcmc", help="The total number of steps in the Markov chain. Default: 1000", default=False);
    parser.add_argument("-chain", dest="chain", help="The number of chains. Default: 1", default=False);
    # MCMC options

    parser.add_argument("-path", dest="phyloacc_path", help="The path to the PhyloAcc binary. Default: PhyloAcc", default=False);
    # Dependency paths
    ## Note: For now we will likely need three dependency paths for the models, but eventually these should all be consolidated
    
    parser.add_argument("-n", dest="num_procs", help="The number of processes that this script should use. Default: 1.", type=int, default=1);
    parser.add_argument("-p", dest="phyloacc_procs", help="The total number of processes that PhyloAcc can spawn. Should be set to the number of expected threads available. Default: 1.", type=int, default=1);
    parser.add_argument("-j", dest="num_jobs", help="The number of jobs (batches) to run in parallel. Must be less than or equal to the total processes for PhyloAcc (-p). Default: 1.", type=int, default=1);
    # User params

    parser.add_argument("-batch", dest="batch_size", help="The number of loci to run per batch. Default: 50", default=False);
    # Batch options

    parser.add_argument("-part", dest="cluster_part", help="The partition or list of partitions (separated by commas) on which to run PhyloAcc jobs.", default=False);
    parser.add_argument("-nodes", dest="cluster_nodes", help="The number of nodes on the specified partition to submit jobs to. Default: 1.", default=False);
    parser.add_argument("-mem", dest="cluster_mem", help="The max memory for each job in GB. Default: 4.", default=False);
    parser.add_argument("-time", dest="cluster_time", help="The time in hours to give each job. Default: 1.", default=False);
    # Cluster options
    
    parser.add_argument("--labeltree", dest="labeltree", help="Simply reads the tree from the input mod file (-m), labels the internal nodes, and exits.", action="store_true", default=False);
    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    # User options
    
    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    parser.add_argument("--depcheck", dest="depcheck", help="Run this to check that all dependencies are installed at the provided path. No other options necessary.", action="store_true", default=False);
    #parser.add_argument("--dryrun", dest="dryrun", help="With all options provided, set this to run through the whole pseudo-it pipeline without executing external commands.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent PhyloAcc from reporting detailed information about each step.", action="store_true", default=False);
    # Run options
    
    parser.add_argument("--qstats", dest="qstats", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    
    args = parser.parse_args();
    # The input options and help messages

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    if args.info_flag:
        globs['info'] = True;
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    globs['label-tree'] = args.labeltree;
    # Parse the --labeltree option

    if not globs['label-tree']:
        globs, deps_passed = PC.execCheck(globs, args);
        if args.depcheck:
            if deps_passed:
                print("\n# All dependencies PASSED.\n")
                sys.exit(0);
            else:
                print("\n# Some dependencies NOT FOUND. Please check your installations and provided paths.\n");
                sys.exit(1);
        # Check the dependency paths

        if args.norun:
            globs['norun'] = True;
            globs['log-v'] = -1;
        globs['overwrite'] = args.ow_flag;
        # Check run mode options.

        if (not args.aln_file and not args.aln_dir) or (args.aln_file and args.aln_dir):
            PC.errorOut("OP1", "One input method must be specified: -a or -d", globs);
        # Check that only one input type is specified

        if args.aln_file:
            if not args.bed_file:
                PC.errorOut("OP2", "A bed file with locus coordinates must be specified with -b when a concatenated alignment file is given with -a", globs);
            globs['aln-file'] = args.aln_file;
            globs['bed-file'] = args.bed_file;
            globs['id-file'] = args.id_file;
        elif args.aln_dir:
            globs['aln-dir'] = args.aln_dir;
        # Save the input type as a global param

    if not args.mod_file:
        PC.errorOut("OP4", "A mod file must be provided with -m", globs);
    globs['mod-file'] = args.mod_file;
    # Check the mod file

    globs = PC.fileCheck(globs);
    # Make sure all the input file actually exist, and get their
    # full paths

    for line in open(globs['mod-file']):
        if line.startswith("TREE: "):
            globs['tree-string'] = line.strip().replace("TREE: ", "");
    # Read the tree string from the MOD file

    try:
        globs['tree-dict'], globs['labeled-tree'], globs['root-node'] = TREE.treeParse(globs['tree-string']);
        globs['tree-tips'] = [ n for n in globs['tree-dict'] if globs['tree-dict'][n][2] == "tip" ];
    except:
        PC.errorOut("OP5", "Error reading tree from mod file!", globs);
    # Read the tree as a dictionary

    if args.labeltree:
        print("# --labeltree SET. LABELING INPUT TREE AND EXITING:\n")
        print(globs['labeled-tree']);
        sys.exit(0);
    # If --labeltree is set, print the tree here and exit

    if not args.targets:
        PC.errorOut("OP6", "Target (-t) species must be specified.", globs);       
    globs['targets'] = args.targets.replace("; ", ";").split(";");
    if args.outgroup:
        globs['outgroup'] = args.outgroup.replace("; ", ";").split(";");
    # Read the target and outgroup groups (if provided)

    if args.conserved:
        globs['conserved'] = args.conserved.replace("; ", ";").split(";");
    else:
        globs['conserved'] = [ node for node in globs['tree-dict'] 
                                if globs['tree-dict'][node][2] =='tip' 
                                and node not in globs['targets'] 
                                and node not in globs['outgroup'] ];
    # Read the conserved group if provided, and if not infer it from the other groups

    for group in ['targets', 'outgroup', 'conserved']:
        for species in globs[group]:
            if species not in globs['tree-dict']:
                PC.errorOut("OP7", "The following species label was provided in a group but is not present in the tree: " + species, globs);
    # A preliminary check here to make sure all provided labels are actually in the input tree

    opt_keys = {'burnin' : args.burnin, 'mcmc' : args.mcmc, 'chain' : args.chain };
    for opt in opt_keys:
        if opt_keys[opt]:
            if not PC.isPosInt(opt_keys[opt]):
                PC.errorOut("OP8", "-b, -m, and -n must all be positive integers.", globs);
            globs[opt] = opt_keys[opt];
    # Get the MCMC options

    if not args.out_dest:
        globs['outdir'] = "phyloacc-out-" + globs['startdatetime'];
    else:
        globs['outdir'] = args.out_dest;

    if not globs['overwrite'] and os.path.exists(globs['outdir']):
        PC.errorOut("OP9", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);

    if not os.path.isdir(globs['outdir']) and not globs['norun']:
        os.makedirs(globs['outdir']);
    # Main output dir

    globs['job-dir'] = os.path.join(globs['outdir'], "phyloacc-job-files");
    if not os.path.isdir(globs['job-dir']) and not globs['norun']:
        os.makedirs(globs['job-dir']);
    # Main job file dir

    job_sub_dirs = { 'job-alns' : 'alns', 'job-cfgs' : 'cfgs', 'job-bed' : 'bed', 'job-smk' : 'snakemake', 'job-out' : 'phyloacc-output' };
    if not globs['norun']:
        for subdir in job_sub_dirs:
            globs[subdir] = os.path.join(globs['job-dir'], job_sub_dirs[subdir]);
            if not os.path.isdir(globs[subdir]):
                os.makedirs(globs[subdir]);
    # Job output subdirs    

    globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
    globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
    # Log file

    if args.batch_size:
        print(args.batch_size);
        if not PC.isPosInt(args.batch_size):
            PC.errorOut("OP10", "The number of loci per batch (-batch) must be a positive integer.", globs);
        else:
            globs['batch-size'] = int(args.batch_size);
    # Batch size

    globs['phyloacc-procs'] = PC.isPosInt(args.phyloacc_procs, default=1);
    globs['num-jobs'] = PC.isPosInt(args.num_jobs, default=1);
    if globs['num-jobs'] > globs['phyloacc-procs']:
        PC.errorOut("OP11", "The specified number of jobs (-j) should not exceed the specified number of processes (-p).", globs);
    globs['procs-per-job'] = math.floor(globs['phyloacc-procs'] / globs['num-jobs'])
    # Determine resource allocation for PhyloAcc

    if not args.cluster_part:
        PC.errorOut("OP12", "At least one cluster partition must be specified with -part.", globs);
    else:
        globs['partition'] = args.cluster_part;
    # Cluster partition option (required)

    if args.cluster_nodes:
        if not PC.isPosInt(args.cluster_nodes):
            PC.errorOut("OP13", "The number of nodes specified (-nodes) must be a positive integer.", globs);
        else:
            globs['num-nodes'] = args.cluster_nodes;
    # Cluster memory option

    if args.cluster_mem:
        if not PC.isPosInt(args.cluster_mem):
            PC.errorOut("OP14", "The specified cluster memory (-mem) must be a positive integer in GB.", globs);
        else:
            globs['mem'] = args.cluster_mem;
    # Cluster memory option

    if args.cluster_time:
        if not PC.isPosInt(args.cluster_time):
            PC.errorOut("OP15", "The specified cluster time (-time) must be a positive integer in hours.", globs);
        else:
            globs['time'] = args.cluster_time + ":00:00";
    # Cluster memory option

    if args.quiet_flag:
        globs['quiet'] = True;
    # Check the quiet option

    globs['num-procs'] = PC.isPosInt(args.num_procs, default=1);
    globs['aln-pool'] = mp.Pool(processes=globs['num-procs']);
    globs['scf-pool'] = mp.Pool(processes=globs['num-procs']);
    # Create the pool of processes for sCF calculation here so we copy the memory profile of the parent process
    # before we've read any large data in

    if args.qstats:
        globs['qstats'] = True;
    # Check for the internal quartet stats option to write to a file.

    if globs['psutil']:
        globs['pids'] = [psutil.Process(os.getpid())];
    # Get the starting process ids to calculate memory usage throughout.

    startProg(globs);
    # After all the essential options have been set, call the welcome function.

    return globs;

#############################################################################

def startProg(globs):
# A nice way to start the program.
    print("#");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Welcome to PhyloAcc -- Bayesian rate analysis of conserved non-coding genomic elements.");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Version (interface) " + globs['interface-version'] + " released on " + globs['releasedate']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# PhyloAcc was developed by Zhirui Hu, Han Yan, Taehee Lee, Gregg Thomas, Tim Sackton, Scott Edwards, and Jun Liu");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Citation:      " + globs['doi']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Report issues: " + globs['github']);
    PC.printWrite(globs['logfilename'], globs['log-v'], "#");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is: " + PC.getDateTime());
    PC.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:              " + globs['pyver'] + "\n#");
    PC.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 45;
    opt_pad = 50;
    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# INPUT/OUTPUT INFO:");

    if globs['aln-file']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Alignment file:", pad) + globs['aln-file']);
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Bed file:", pad) + globs['bed-file']);
    elif globs['aln-dir']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Alignment directory:", pad) + globs['aln-dir']);

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Tree/rate file (mod file from PHAST):", pad) + globs['mod-file']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Tree read from mod file:", pad) + globs['tree-string']);

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Output directory:", pad) + globs['outdir']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc run directory:", pad) + globs['job-dir']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));
    # Input/Output
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# DEPENDENCY PATHS:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Program", pad) + "Specified Path");

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc", pad) + globs['phyloacc']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# PhyloAcc-gBGC", pad) + globs['phyloacc-gbgc']);
    # Dependency paths
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# SPECIES GROUPS:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Group", pad) + PC.spacedOut("Species", opt_pad));

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Targets (-t)", pad) + ";".join(globs['targets'])); 
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Conserved (-c)", pad) + ";".join(globs['conserved']));
    if globs['outgroup']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Outgroups (-g)", pad) + ";".join(globs['outgroup']));
    # Species groups
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# CLUSTER OPTIONS:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Option", pad) + PC.spacedOut("Setting", opt_pad));

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Partition(s)", pad) + globs['partition']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Number of nodes", pad) + globs['num-nodes']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Max mem per job (gb)", pad) + globs['mem']);
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Time per job", pad) + globs['time']); 
    # Cluster options
    #######################

    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Option", pad) + PC.spacedOut("Current setting", opt_pad) + "Current action");

    if globs['aln-file']:
        if globs['id-file']:
            PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -i:", pad) + 
                PC.spacedOut(str(globs['id-file']), opt_pad) +
                "Only loci names specified in this file will be tested.");
        else:
            PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -i:", pad) + 
                PC.spacedOut(str(globs['id-file']), opt_pad) +
                "No ID file provided, all loci in input bed file will be tested.");           

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -burnin:", pad) + 
                PC.spacedOut(str(globs['burnin']), opt_pad) +
                "This number of steps in the chain will discarded as burnin");         
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -mcmc:", pad) + 
                PC.spacedOut(str(globs['mcmc']), opt_pad) +
                "The number of steps in each chain");      
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -chain:", pad) + 
                PC.spacedOut(str(globs['chain']), opt_pad) +
                "The number of chains to run");      
    # MCMC options

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Loci per batch (-batch)", pad) + 
                PC.spacedOut(str(globs['batch-size']), opt_pad) + 
                "PhyloAcc will run this many loci in a single command.");
    # Batch size

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Processes (-p)", pad) + 
                PC.spacedOut(str(globs['num-procs']), opt_pad) + 
                "PhyloAcc and this interface will use this many total processes.");
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Jobs (-j)", pad) + 
                PC.spacedOut(str(globs['num-jobs']), opt_pad) + 
                "PhyloAcc will spawn this many jobs in parallel.");
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Processes per job", pad) + 
                PC.spacedOut(str(globs['procs-per-job']), opt_pad) + 
                "Each job will use this many processes.");

    if globs['num-procs'] % globs['num-jobs'] != 0:
        PC.printWrite(globs['logfilename'], globs['log-v'], 
            "# WARNING: The number of jobs is not a multiple of the number of processes! Some processes will be idle.");
    # Reporting the resource options.

    if globs['overwrite']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --overwrite", pad) +
                    PC.spacedOut("True", opt_pad) + 
                    "PhyloAcc will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option.

    if not globs['quiet']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --quiet", pad) + 
                    PC.spacedOut("False", opt_pad) + 
                    "Time, memory, and status info will be printed to the screen while PhyloAcc is running.");
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --quiet", pad) + 
                    PC.spacedOut("True", opt_pad) + 
                    "No further information will be printed to the screen while PhyloAcc is running.");
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
        PC.printWrite(globs['logfilename'], globs['log-v'], "# Running...");
    # Reporting the quiet option.

    # if globs['debug']:
    #     PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --debug", pad) + 
    #                 PC.spacedOut("True", opt_pad) + 
    #                 "Printing out a bit of debug info.");
    # Reporting the debug option.

    if globs['qstats']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --qstats", pad) + 
            PC.spacedOut("True", opt_pad) + 
            "Writing out a file with quartet site counts.");

    if globs['norun']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --norun", pad) + 
                    PC.spacedOut("True", opt_pad) + 
                    "ONLY PRINTING RUNTIME INFO.");
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    # Reporting the norun option.

    # Other options
    #######################

#############################################################################