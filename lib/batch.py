#############################################################################
# Functions to generate files for batch jobs of PhyloAcc (with snakemake)
# Gregg Thomas
#############################################################################

import os
import lib.core as PC

#############################################################################

def genJobFiles(globs):
# This function generates all files for phyloacc to run on a per locus basis with snakemake

    step = "Writing PhyloAcc job files";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    jobs = 0;
    for aln in globs['alns']:
    # For every input alignment, write an alignment file, config, file, bed, file and
    # output directory.
    ## TODO: Parallelize?
    ## TODO: Don't need to re-write individual aligns if that is how input is provided...

        cur_out_dir = os.path.join(globs['job-out'], aln + "-phyloacc-out");
        if not os.path.isdir(cur_out_dir):
            os.makedirs(cur_out_dir);
        # Make the phyloacc output directory for the current alignment

        aln_len = globs['aln-stats'][aln]['length']
        # Get the alignment length for the bed file

        cur_aln_file = os.path.join(globs['job-alns'], aln + ".fa");
        with open(cur_aln_file, "w") as alnfile:
            for spec in globs['alns'][aln]:
                alnfile.write(">" + spec + "\n");
                alnfile.write(globs['alns'][aln][spec] + "\n");
        # Write the current alignment to file
        
        cur_bed_file = os.path.join(globs['job-bed'], aln + ".bed");
        with open(cur_bed_file, "w") as bedfile:
            outline = [aln, "0", str(aln_len), aln];
            bedfile.write("\t".join(outline));   
        # Write a bed file that contains all coordinates for the current alignment
        # Only need to do this because the PhyloAcc config file requires a SEG_FILE
        
        cur_cfg_file = os.path.join(globs['job-cfgs'], aln + ".cfg");
        with open(cur_cfg_file, "w") as cfgfile:
            cfgfile.write("PHYTREE_FILE " + os.path.abspath(globs['mod-file']) + "\n");
            cfgfile.write("SEG_FILE " + os.path.abspath(cur_bed_file) + "\n");
            cfgfile.write("ALIGN_FILE " + os.path.abspath(cur_aln_file) + "\n");
            cfgfile.write("RESULT_FOLDER " + os.path.abspath(cur_out_dir) + "\n");
            cfgfile.write("PREFIX " + aln + "\n");
            cfgfile.write("BURNIN " + str(globs['burnin']) + "\n");
            cfgfile.write("MCMC " + str(globs['mcmc']) + "\n");
            cfgfile.write("CHAIN " + str(globs['chain']) + "\n");
            cfgfile.write("TARGETSPECIES " + str(globs['targets']) + "\n");
            cfgfile.write("OUTGROUP " + str(globs['outgroup']) + "\n");
            cfgfile.write("CONSERVE " + str(globs['conserved']) + "\n");
            cfgfile.write("GAPCHAR -\n");
            cfgfile.write("NUM_THREAD " + str(globs['procs-per-job']) + "\n");
            cfgfile.write("VERBOSE 1");
        # Write the phyloacc config file for the current alignment

        jobs += 1;
        # Iterate the number of jobs successfully written for the status update

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(jobs) + " jobs written");
    # Status update       

    return(globs);

#############################################################################

def writeSnakemake(globs):

    step = "Writing Snakemake config file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    globs['smk-config'] = os.path.join(globs['job-smk'], "phyloacc-config.yaml");

    with open(globs['smk-config'], "w") as configfile:
        configfile.write("input_directory: " + os.path.abspath(globs['job-cfgs']) + "\n");
        configfile.write("output_directory: " + os.path.abspath(globs['job-out']) + "\n");
        configfile.write("locus_list: " + str(list(globs['alns'].keys())) + "\n");

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake config written");
    # Status update  

    step = "Writing Snakemake cluster profile";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    globs['profile-dir'] = os.path.join(globs['job-smk'], "profiles", "slurm_profile");
    if not os.path.isdir(globs['profile-dir']):
        os.makedirs(globs['profile-dir']);

    profile_file = os.path.join(globs['profile-dir'], "config.yaml");
    ## This is a profile for SLURM. Will likely need templates for different job schedulers

    cluster_logdir = os.path.abspath(os.path.join(globs['job-dir'], "slurm-logs"));

    with open(profile_file, "w") as pro_file:
        pro_file.write("jobs: " + str(globs['num-jobs']) + "\n");
        pro_file.write("cluster:\n");
        pro_file.write("  mkdir -p " + cluster_logdir + "/{rule}/ &&\n");
        pro_file.write("  sbatch\n");
        pro_file.write("  --partition={resources.partition}\n");
        pro_file.write("  --nodes={resources.nodes}\n");
        pro_file.write("  --cpus-per-task=" + str(globs['procs-per-job']) + "\n");
        pro_file.write("  --job-name={rule}-{wildcards}\n");
        pro_file.write("  --mem={resources.mem}\n");
        pro_file.write("  --time={resources.time}\n");
        pro_file.write("  --output=" + cluster_logdir + "/{rule}/{rule}-{wildcards}-%j.out\n");
        #pro_file.write("  --mail-type=END,FAIL\n");
        pro_file.write("  --mail-user=gthomas@g.harvard.edu\n");
        pro_file.write("default-resources:\n");
        pro_file.write("  - partition='" + globs['partition'] + "'\n");
        pro_file.write("  - nodes='" + globs['num-nodes'] + "'\n");
        pro_file.write("  - mem='" + globs['mem'] + "g'\n");
        pro_file.write("  - time='" + globs['time'] + "'\n");
        pro_file.write("latency-wait: 30\n");
        pro_file.write("verbose: true");

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake profile written");
    # Status update     

    return globs; 

#############################################################################