#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import argparse
import lib.core as PC

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py
    parser = argparse.ArgumentParser(description="PhyloAcc: Bayesian rate analysis of conserved non-coding genomic elements");

    parser.add_argument("-a", dest="aln_file", help="An alignment file with all loci concatenated. -b and -i must also be specified. Expected as FASTA format for now.", default=False);
    parser.add_argument("-b", dest="bed_file", help="A bed file with coordinates for the loci in the concatenated alignment file. -a and -i must also be specified.", default=False);
    parser.add_argument("-i", dest="id_file", help="A text file with locus names, one per line, corresponding to regions in the input bed file. -a and -b must also be specified.", default=False);

    parser.add_argument("-d", dest="aln_dir", help="A directory containing individual alignment files for each locus. Expected as FASTA format for now.", default=False);
    
    parser.add_argument("-m", dest="mod_file", help="A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from PHAST.", default=False);
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: phyloacc-[date]-[time]", default=False);
    # Output
    
    parser.add_argument("-t", dest="targets", help="Tip labels in the input tree to be used as target species. Enter multiple labels separated by semi-colons (;).", default=False);
    parser.add_argument("-c", dest="conserved", help="Tip labels in the input tree to be used as conserved species. Enter multiple labels separated by semi-colons (;).", default=False);
    parser.add_argument("-g", dest="outgroup", help="Tip labels in the input tree to be used as outgroup species. Enter multiple labels separated by semi-colons (;).", default=False);
    # Phylo options
    
    parser.add_argument("-burnin", dest="burnin", help="The number of steps to be discarded in the Markov chain as burnin. Default: 500", default=False);
    parser.add_argument("-mcmc", dest="mcmc", help="The total number of steps in the Markov chain. Default: 1000", default=False);
    parser.add_argument("-chain", dest="chain", help="The number of chains. Default: 1", default=False);
    # MCMC options

    parser.add_argument("-p", dest="phyloacc_path", help="The path to the PhyloAcc binary. Default: PhyloAcc", default=False);
    # Dependency paths
    ## Note: For now we wil likely need three dependency paths for the models, but eventually these should all be consolidated
    parser.add_argument("-j", dest="num_jobs", help="The number of jobs (loci) to run in parallel. Should be set to the number of expected threads available. Default: 1.", type=int, default=1);
    # User params
    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    # User options
    parser.add_argument("--depcheck", dest="depcheck", help="Run this to check that all dependencies are installed at the provided path. No other options necessary.", action="store_true", default=False);
    #parser.add_argument("--dryrun", dest="dryrun", help="With all options provided, set this to run through the whole pseudo-it pipeline without executing external commands.", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent PhyloAcc from reporting detailed information about each step.", action="store_true", default=False);
    # Run options
    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    args = parser.parse_args();
    # The input options and help messages

    globs, deps_passed = PC.execCheck(globs, args);
    if args.depcheck:
        if deps_passed:
            print("\n# All dependencies PASSED.\n")
            sys.exit(0);
        else:
            print("\n# Some dependencies NOT FOUND. Please check your installations and provided paths.\n");
            sys.exit(1);
    # Check the dependency paths.

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
        if not args.id_file:
            PC.errorOut("OP3", "An ID file with locus names must be specified with -i when a concatenated alignment file is given with -a", globs);
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

    if not all(labels for labels in [args.targets, args.conserved, args.outgroup]):
        PC.errorOut("OP5", "Target (-t), conserved (-c), and outgroup (-g) species must all be specified.", globs);
    globs['targets'] = args.targets.replace("; ", ";").split(";");
    globs['conserved'] = args.conserved.replace("; ", ";").split(";");
    globs['outgroup'] = args.outgroup.replace("; ", ";").split(";");
    # Read the species groups
    ## May also want to do a preliminary check here to make sure all provided labels are actually in the input tree
    ## Would also require reading the tree somewhere around here

    opt_keys = {'burnin' : args.burnin, 'mcmc' : args.mcmc, 'chain' : args.chain };
    for opt in opt_keys:
        if opt_keys[opt]:
            if not PC.isPosInt(opt_keys[opt]):
                PC.errorOut("OP6", "-b, -m, and -n must all be positive integers.", globs);
            globs[opt] = opt_keys[opt];
    # Get the MCMC options

    if not args.out_dest:
        globs['outdir'] = "phyloacc-out-" + globs['startdatetime'];
    else:
        globs['outdir'] = args.out_dest;

    if not globs['overwrite'] and os.path.exists(globs['outdir']):
        PC.errorOut("OP7", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);

    if not os.path.isdir(globs['outdir']) and not globs['norun']:
        os.makedirs(globs['outdir']);
    # Main output dir

    globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
    globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
    # Log file

    globs['jobs'] = PC.isPosInt(args.num_jobs);
    # Number of jobs

    if args.quiet_flag:
        globs['quiet'] = True;
    # Check the quiet option

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
    PC.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + " ".join(sys.argv) + "\n#");

    #######################

    pad = 45;
    opt_pad = 50;
    PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    PC.printWrite(globs['logfilename'], globs['log-v'], "# INPUT/OUTPUT INFO:");

    if globs['aln-file']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Alignment file:", pad) + globs['aln-file']);
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Bed file:", pad) + globs['bed-file']);
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# ID file:", pad) + globs['id-file']);
    elif globs['aln-dir']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Alignment directory:", pad) + globs['aln-dir']);

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Tree/rate file (mod file from PHAST):", pad) + globs['mod-file']);

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Output directory:", pad) + globs['outdir']);
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
    PC.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");    
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Option", pad) + PC.spacedOut("Current setting", opt_pad) + "Current action");

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -t", pad) +
                PC.spacedOut(";".join(globs['targets']), opt_pad) +
                "These are the target species"); 
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -c", pad) + 
                PC.spacedOut(";".join(globs['conserved']), opt_pad) +
                "These are the conserved species");
    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -g:", pad) + 
                PC.spacedOut(";".join(globs['outgroup']), opt_pad) +
                "These are the outgroup species");
    # Species groups

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

    if globs['overwrite']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --overwrite", pad) +
                    PC.spacedOut("True", opt_pad) + 
                    "PhyloAcc will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option.

    PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -j", pad) + 
                PC.spacedOut(str(globs['num-jobs']), opt_pad) + 
                "PhyloAcc will spawn this many jobs in parallel.");
    # Reporting the processes option.

    if not globs['quiet']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --quiet", pad) + 
                    PC.spacedOut("False", opt_pad) + 
                    "Runtime, memory, and command info will be printed to the screen while PhyloAcc is running.");
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

    if globs['norun']:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --norun", pad) + 
                    PC.spacedOut("True", opt_pad) + 
                    "ONLY PRINTING RUNTIME INFO.");
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    # Reporting the norun option.

    # Other options
    #######################