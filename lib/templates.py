#############################################################################
# Templates for various files written by the phyloacc interface
#############################################################################

def phyloaccConfig():

    phyloacc_template = """PHYTREE_FILE {mod_file}
SEG_FILE {bed_file}
ALIGN_FILE {aln_file}
RESULT_FOLDER {outdir}
PREFIX {batch}
BURNIN {burnin}
MCMC {mcmc}
CHAIN {chain}
TARGETSPECIES {targets}
OUTGROUP {outgroup}
CONSERVE {conserved}
GAPCHAR -
NUM_THREAD {procs_per_job}
VERBOSE 1
"""

    return phyloacc_template;

#############################################################################

def snakemake():

    smk_template = """#############################################################################
# Pipeline for running phyloacc per locus
# Generated from: {cmd}
# On: {dt}
#############################################################################

#############################################################################

import os 

#############################################################################

INDIR = config["input_directory"];
OUTDIR = config["output_directory"];
BATCHES = config["batch_list"];
# Inputs for the snakemake pipeline are read from the config file generated by
# the interface

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(OUTDIR, "{{batch}}-phyloacc-out", "{{batch}}_elem_lik.txt"), batch=BATCHES)
# This rule just checks that some of the phyloacc output files are present for each batch

#############################################################################
# Pipeline rules

rule run_phyloacc:
    input:
        os.path.join(INDIR, "{{batch}}.cfg")
    output:
        os.path.join(OUTDIR, "{{batch}}-phyloacc-out", "{{batch}}_elem_lik.txt")
    log:
        os.path.join(OUTDIR, "{{batch}}-phyloacc-out", "{{batch}}-phyloacc.log")
    shell:
        \"\"\"
        {path} {{input}} &> {{log}}
        \"\"\"
        # Need to replace this with path to PhyloAcc on current install... depends on how we package the codebases together
        
# This rule runs phyloacc on each batch individually. With a cluster profile
# each batch will be submitted as a job.

#############################################################################
"""

    return smk_template;

#############################################################################

def snakemakeConfig():
    
    config_template = """input_directory: {indir}
output_directory: {outdir}
batch_list: {batches}
"""

    return config_template;

#############################################################################

def snakemakeProfile():
    
    profile_template = """jobs: {num_jobs}
cluster:
  mkdir -p {cluster_logdir}/{{rule}}/ &&
  sbatch
  --partition={{resources.partition}}
  --nodes={{resources.nodes}}
  --cpus-per-task={procs_per_job}
  --job-name={{rule}}-{{wildcards}}
  --mem={{resources.mem}}
  --time={{resources.time}}
  --output={cluster_logdir}/{{rule}}/{{rule}}-{{wildcards}}-%j.out
default-resources:
  - partition='{part}'
  - nodes='{num_nodes}'
  - mem='{mem}g'
  - time='{time}'
latency-wait: 30
verbose: true
"""

    return profile_template;

#############################################################################