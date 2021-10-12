#############################################################################
# Functions to generate files for batch jobs of PhyloAcc (with snakemake)
# Gregg Thomas
#############################################################################

import os
import lib.core as PC
import lib.templates as TEMPLATES

#############################################################################

def genJobFiles(globs):
# This function generates all files for phyloacc to run on a per locus basis with snakemake

    step = "Writing PhyloAcc job files";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    batch_num = 0;
    for batch in PC.dictChunk(globs['alns'], globs['batch-size']):
    # Split the sequence dictionary by batches and yield the result here

        batch_num += 1;
        batch_num_str = str(batch_num);
        # Batch counting

        cur_out_dir = os.path.join(globs['job-out'], batch_num_str + "-phyloacc-out");
        if not os.path.isdir(cur_out_dir):
            os.makedirs(cur_out_dir);
        # Make the phyloacc output directory for the current batch

        batch_concat = { label : "" for label in globs['tree-tips'] };
        # Dictionary to concatenate alignments for the current batch

        batch_aln_list = [];
        # As we concatenate, add the alignment name to this list so we can go through
        # in the same order as we make the bed file

        for aln in batch:
            #print(aln);
            batch_aln_list.append(aln);
            for label in globs['tree-tips']:
                batch_concat[label] += batch[aln][label];
            # May need a check here for missing sequences
        # Concatenate alignments in the current batch together

        cur_aln_file = os.path.join(globs['job-alns'], batch_num_str + ".fa");
        with open(cur_aln_file, "w") as alnfile:
            for spec in globs['tree-tips']:
                alnfile.write(">" + spec + "\n");
                alnfile.write(batch_concat[spec] + "\n");
        # Write the current concatenated alignment to file

        len_sum = 0;
        # Keeps track of the last coordinate written in the bed file

        cur_bed_file = os.path.join(globs['job-bed'], batch_num_str + ".bed");
        with open(cur_bed_file, "w") as bedfile:
            for aln in batch_aln_list:

                aln_len = globs['aln-stats'][aln]['length'];
                # Get the alignment length for the bed file

                end_coord = len_sum + aln_len;
                # Get the end coordinate of the current locus in the concatenated alignment by
                # adding the length to the previous length sum

                outline = [aln, str(len_sum), str(end_coord), aln];
                bedfile.write("\t".join(outline) + "\n");
                # Write the info for the current locus

                len_sum += aln_len;
                # Add the length to the length sum as the start coordinate for the next locus
        # Write a bed file that contains all coordinates for the current alignment

        cur_cfg_file = os.path.join(globs['job-cfgs'], batch_num_str + ".cfg");
        with open(cur_cfg_file, "w") as cfgfile:
            cfgfile.write(TEMPLATES.phyloaccConfig().format(mod_file=os.path.abspath(globs['mod-file']),
                                                                    bed_file=os.path.abspath(cur_bed_file),
                                                                    aln_file=os.path.abspath(cur_aln_file),
                                                                    outdir=os.path.abspath(cur_out_dir),
                                                                    batch=batch_num_str,
                                                                    burnin=str(globs['burnin']),
                                                                    mcmc=str(globs['mcmc']),
                                                                    chain=str(globs['chain']),
                                                                    targets= ";".join(globs['targets']),
                                                                    outgroup=";".join(globs['outgroup']),
                                                                    conserved=";".join(globs['conserved']),
                                                                    procs_per_job=str(globs['procs-per-job'])
                                                                    ))
        # Write the phyloacc config file for the current concatenated alignment
    ## End batch loop

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(batch_num) + " jobs written");
    # Status update

    globs['num-batches'] = batch_num;
    # Can only compute this here after we've split the number of input alignments by batch size

    return(globs);

############################################################################# 

def writeSnakemake(globs):
# A function to write the various snakemake files

    step = "Writing Snakemake file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update
    
    globs['smk'] = os.path.join(globs['job-smk'], "run_phyloacc.smk");

    with open(globs['smk'], "w") as smkfile:
        smkfile.write(TEMPLATES.snakemake().format(cmd=globs['call'],
                                                       dt=PC.getDateTime(),
                                                       path=os.path.abspath(globs['phyloacc'])
                                                       ))

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake file written");
    # Status update  

    ####################

    step = "Writing Snakemake config file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    globs['smk-config'] = os.path.join(globs['job-smk'], "phyloacc-config.yaml");

    with open(globs['smk-config'], "w") as configfile:
        configfile.write(TEMPLATES.snakemakeConfig().format(indir=os.path.abspath(globs['job-cfgs']),
                                                            outdir=os.path.abspath(globs['job-out']),
                                                            batches=str(list(range(1, globs['num-batches']+1)))
                                                            ))

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake config written");
    # Status update  

    ####################

    step = "Writing Snakemake cluster profile";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status update

    globs['profile-dir'] = os.path.join(globs['job-smk'], "profiles", "slurm_profile");
    if not os.path.isdir(globs['profile-dir']):
        os.makedirs(globs['profile-dir']);

    profile_file = os.path.join(globs['profile-dir'], "config.yaml");
    ## This is a profile for SLURM. Will likely need templates for different job schedulers

    cluster_logdir = os.path.abspath(os.path.join(globs['job-dir'], "slurm-logs"));

    with open(profile_file, "w") as profile:
        profile.write(TEMPLATES.snakemakeProfile().format(num_jobs=str(globs['num-jobs']),
                                                          cluster_logdir=cluster_logdir,
                                                          procs_per_job=str(globs['procs-per-job']),
                                                          part=globs['partition'],
                                                          num_nodes=globs['num-nodes'],
                                                          mem=globs['mem'],
                                                          time=globs['time']
                                                        ))

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: Snakemake profile written");
    # Status update     

    return globs; 

#############################################################################