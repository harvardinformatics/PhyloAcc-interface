#############################################################################
# Functions to generate plots and html files for easy visualization of 
# input dataset
# Gregg Thomas
#############################################################################

import os
import shutil
import re
import lib.core as PC
import lib.tree as TREE
import lib.templates_2 as TEMPLATES
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D as lines
import matplotlib.patches as mpatches

from Bio import Phylo
from io import StringIO

#############################################################################

def genPlots(globs):

    step = "Generating summary plots";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status updated

    mpl.rcParams["axes.spines.right"] = False;
    mpl.rcParams["axes.spines.top"] = False;

    mpl.rcParams["axes.labelsize"] = 16;
    mpl.rcParams["axes.edgecolor"] = "#595959";
    mpl.rcParams["axes.linewidth"] = 1.5;

    mpl.rcParams["xtick.labelcolor"] = "#595959";
    mpl.rcParams["xtick.labelsize"] = 14;

    mpl.rcParams['xtick.color'] = "#595959";
    mpl.rcParams['xtick.major.size'] = 6;
    mpl.rcParams['xtick.major.width'] = 1.5;

    mpl.rcParams["ytick.labelcolor"] = "#595959";
    mpl.rcParams["ytick.labelsize"] = 14;

    mpl.rcParams['ytick.color'] = "#595959";
    mpl.rcParams['ytick.major.size'] = 6;
    mpl.rcParams['ytick.major.width'] = 1.5;
    # Global theme settings for all matplotlib figures following
    
    ####################

    st_file = os.path.join(globs['plot-dir'], globs['input-tree-plot-file']);
    # The file to save the species tree figure

    branch_cols = PC.coreCol(pal="wilke", numcol=3);
    # The colors for the target, conserved, and outgroup branches

    num_spec = len(globs['tree-tips']);
    # The number of species in the input free, to adjust height of figure

    tree_str = TREE.addBranchLength(globs['labeled-tree'], globs['tree-dict'], no_label=True);
    #print(tree_str);
    # Re-add branch lengths and remove labels to the input tree for plotting

    handle = StringIO(tree_str);
    tree = Phylo.read(handle, "newick");
    # Parse the tree string with Bio
    
    for clade in tree.get_terminals():
        if clade.name in globs['targets']:
            clade.color = branch_cols[0];
            target_lg = lines([0], [0], label="Targets", color = branch_cols[0]);
        elif clade.name in globs['conserved']:
            clade.color = branch_cols[1];
            conserve_lg = lines([0], [0], label="Conserved", color = branch_cols[1]);
        elif clade.name in globs['outgroup']:
            clade.color = branch_cols[2];
            outgroup_lg = lines([0], [0], label="Outgroup", color = branch_cols[2]);
    # Color the tip branches based on their input category and specify their legend entries

    fig = plt.figure(figsize=(num_spec/2.54, 25.4/2.54));
    # Specify the plot size depending on the number of species

    axes = fig.add_subplot(1, 1, 1);
    axes.axes.get_yaxis().set_visible(False);
    axes.spines['left'].set_visible(False);
    # Set the axes of the tree figure

    Phylo.draw(tree, axes=axes, show_confidence=False, do_show=False);
    # Draw the tree

    plt.legend(loc='upper left', handles=[target_lg, conserve_lg, outgroup_lg]);
    # Add the legend

    plt.savefig(st_file, dpi=100, bbox_inches='tight');
    # Save the figure

    # Species tree
    ####################

    aln_list = [ aln for aln in globs['aln-stats'] ];
    # A single list of alignment IDs for consistency between dictionary lookups

    aln_len_hist_file = os.path.join(globs['plot-dir'], globs['aln-len-plot-file']);
    aln_lens = [ globs['aln-stats'][aln]['length'] for aln in aln_list ];
    
    plt.figure(figsize=(8,6));
    plt.hist(aln_lens, color=PC.coreCol(pal="wilke", numcol=1, offset=1)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlabel("Alignment length");
    plt.ylabel("# loci");

    plt.savefig(aln_len_hist_file, dpi=100);
    # Locus length (hist)
    ####################

    seq_len_hist_file = os.path.join(globs['plot-dir'], globs['seq-len-plot-file']);
    seq_lens = [ globs['aln-stats'][aln]['avg-nogap-seq-len'] for aln in aln_list ];
    
    plt.figure(figsize=(8,6));
    plt.hist(seq_lens, color=PC.coreCol(pal="wilke", numcol=1, offset=2)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlabel("Avg. sequence length without gaps (bp)");
    plt.ylabel("# loci");

    plt.savefig(seq_len_hist_file, dpi=100);
    # Avg. sequence length without gaps (hist)
    ####################

    # var_sites_hist_file = os.path.join(globs['plot-dir'], "variable-sites-hist.png");
    var_sites = [ globs['aln-stats'][aln]['variable-sites'] for aln in aln_list ];
    
    # plt.figure(figsize=(8,6));
    # plt.hist(var_sites, color=PC.coreCol(pal="wilke", numcol=1)[0], bins="sturges", edgecolor="#999999");
    # #plt.boxplot(locus_len)
    # plt.xlabel("# of variable sites");
    # plt.ylabel("# loci");

    # plt.savefig(var_sites_hist_file);
    # Variable sites (hist)
    ####################

    inf_sites_hist_file = os.path.join(globs['plot-dir'], globs['inf-sites-plot-file']);
    inf_sites = [ globs['aln-stats'][aln]['informative-sites'] for aln in aln_list ];
    
    plt.figure(figsize=(8,6));
    plt.hist(inf_sites, color=PC.coreCol(pal="wilke", numcol=1, offset=3)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlabel("# of informative sites");
    plt.ylabel("# loci");

    plt.savefig(inf_sites_hist_file, dpi=100);
    # Informative sites (hist)
    ####################

    inf_sites_frac_hist_file = os.path.join(globs['plot-dir'], globs['inf-sites-frac-plot-file']);
    #inf_sites_frac = [ globs['aln-stats'][aln]['informative-sites'] / globs['aln-stats'][aln]['length'] for aln in aln_list ];
    inf_sites_frac = [ inf_sites[i] / aln_lens[i] for i in range(len(aln_lens)) ];

    plt.figure(figsize=(8,6));
    plt.hist(inf_sites_frac, color=PC.coreCol(pal="wilke", numcol=1, offset=4)[0], bins="sturges", edgecolor="#999999");
    #plt.boxplot(locus_len)
    plt.xlim([0, 1]);
    plt.xlabel("Fraction of sites that are informative");
    plt.ylabel("# loci");

    plt.savefig(inf_sites_frac_hist_file, dpi=100);
    # Fraction of sites that are informative (hist)
    ####################

    var_inf_sites_file = os.path.join(globs['plot-dir'], globs['var-inf-sites-plot-file']);

    slope, intercept = np.polyfit(var_sites, inf_sites, 1);
    inf_sites_pred = np.polyval([slope, intercept], var_sites);

    plt.figure(figsize=(8,6));
    plt.plot(var_sites, inf_sites_pred, color="#333333", linestyle='dashed', dashes=(5, 20));
    plt.scatter(var_sites, inf_sites, color=PC.coreCol(pal="wilke", numcol=1, offset=5)[0], alpha=0.25);
    plt.xlabel("# of variable sites");
    plt.ylabel("# of informative sites");

    plt.savefig(var_inf_sites_file, dpi=100);
    
    # Variable sites vs. informative sites (scatter w regression)
    ####################

    if globs['run-mode'] == 'adaptive':
        avg_scf_hist_file = os.path.join(globs['plot-dir'], globs['avg-scf-hist-file']);

        avg_scf = [ globs['aln-stats'][aln]['node-scf-avg'] for aln in aln_list if globs['aln-stats'][aln]['node-scf-avg'] != "NA" ];

        plt.figure(figsize=(8,6));
        plt.hist(avg_scf, color=PC.coreCol(pal="wilke", numcol=1, offset=6)[0], bins="sturges", edgecolor="#999999");
        plt.xlim([0, 1]);
        plt.xlabel("Avg. sCF across all branches per locus");
        plt.ylabel("# loci");

        plt.savefig(avg_scf_hist_file, dpi=100);        

        # Avg. scf per locus (hist)
        ####################

        low_scf_hist_file = os.path.join(globs['plot-dir'], globs['low-scf-hist-file']);

        perc_low_scf = [ globs['aln-stats'][aln]['perc-low-scf-nodes'] for aln in aln_list if globs['aln-stats'][aln]['perc-low-scf-nodes'] != "NA" ];

        plt.figure(figsize=(8,6));
        plt.hist(perc_low_scf, color=PC.coreCol(pal="wilke", numcol=1, offset=7)[0], bins="sturges", edgecolor="#999999");
        plt.xlim([0, 1]);
        plt.xlabel("% of branches with sCF below " + str(globs['min-scf']) + " per locus");
        plt.ylabel("# loci");

        plt.savefig(low_scf_hist_file, dpi=100);   

        # % of branches with low sCF per locus (hist)
        ####################

        scf_tree_file = os.path.join(globs['plot-dir'], globs['scf-tree-plot-file']);
        # The file to save the species tree figure

        tree_str = TREE.addBranchLength(globs['labeled-tree'], globs['tree-dict'], no_label=True, keep_tp_label=True);

        for node in globs['scf']:
            tree_str = tree_str.replace(node, node + "_" + str(round(globs['scf'][node]['avg-quartet-scf'], 2)));
        # For every node in the tree, add the averaged scf value over all loci to the label

        tree_str = re.sub("<[\d]+>[_]?", "", tree_str);
        #tree_str = re.sub("<[\d]+>", "", tree_str);

        #print(tree_str);
        # Re-add branch lengths and remove labels to the input tree for plotting

        handle = StringIO(tree_str);
        tree = Phylo.read(handle, "newick");
        # Parse the tree string with Bio
        
        # for clade in tree.get_terminals():
        #     if clade.name in globs['targets']:
        #         clade.color = branch_cols[0];
        #         target_lg = lines([0], [0], label="Targets", color = branch_cols[0]);
        #     elif clade.name in globs['conserved']:
        #         clade.color = branch_cols[1];
        #         conserve_lg = lines([0], [0], label="Conserved", color = branch_cols[1]);
        #     elif clade.name in globs['outgroup']:
        #         clade.color = branch_cols[2];
        #         outgroup_lg = lines([0], [0], label="Outgroup", color = branch_cols[2]);
        # # Color the tip branches based on their input category and specify their legend entries

        fig = plt.figure(figsize=(num_spec/2.54, 25.4/2.54));
        # Specify the plot size depending on the number of species

        axes = fig.add_subplot(1, 1, 1);
        axes.axes.get_yaxis().set_visible(False);
        axes.spines['left'].set_visible(False);
        # Set the axes of the tree figure

        Phylo.draw(tree, axes=axes, do_show=False);
        # Draw the tree

        # plt.legend(loc='upper left', handles=[target_lg, conserve_lg, outgroup_lg]);
        # Add the legend

        plt.savefig(scf_tree_file, dpi=100, bbox_inches='tight');
        # Save the figure        

        #globs['scf-labeled-tree']
        # scf tree (phylo)
        ####################

        bl_scf_file = os.path.join(globs['plot-dir'], globs['bl-scf-plot-file']);

        bls, scfs = [], [];
        for node in globs['tree-dict']:
            if globs['tree-dict'][node][2] != 'internal':
                continue;

            if node not in globs['scf']:
                continue;

            bls.append(float(globs['tree-dict'][node][0]));
            scfs.append(globs['scf'][node]['avg-quartet-scf']);
        # Gets the values out of their tables

        slope, intercept = np.polyfit(bls, scfs, 1);
        scfs_pred = np.polyval([slope, intercept], bls);

        plt.figure(figsize=(8,6));
        plt.plot(bls, scfs_pred, color="#333333", linestyle='dashed', dashes=(5, 20));
        plt.scatter(bls, scfs, color=PC.coreCol(numcol=1)[0], alpha=0.5);
        plt.xlabel("Branch length");
        plt.ylabel("sCF");

        plt.savefig(bl_scf_file, dpi=100);
        
        # Variable sites vs. informative sites (scatter w regression)
        ####################

    step_start_time = PC.report_step(globs, step, step_start_time, "Success");

#############################################################################

def writeHTML(globs):
    step = "Writing HTML summary file";
    step_start_time = PC.report_step(globs, step, False, "In progress...");
    # Status updated

    # if not os.path.isdir(globs['html-dir']):
    #     shutil.copytree( os.path.join(globs['script-dir'], "html"), globs['html-dir'] );
    # Copy the associated html files (stylesheets, images) from the provided template folders

    if globs['run-mode'] == 'adaptive':
        comment_start = "";
        comment_end = "";
    else:
        comment_start = "<!-- This block is only used when run mode (-r) is adaptive";
        comment_end = "-->";
    # Some HTML blocks will be commented out depending on the run mode

    if globs['theta']:
        theta_comment_start = "";
        theta_comment_end = "";
    else:
        theta_comment_start = "<!-- This block is only displayed when --theta is specified";
        theta_comment_end = "-->";        

    if globs['coal-tree-file']:
        coal_tree_comment_start = "";
        coal_tree_comment_end = "";
    else:
        coal_tree_comment_start = "<!-- This block is only displayed when -l is specified";
        coal_tree_comment_end = "-->";            

    with open(globs['html-file'], "w") as htmlfile:
        htmlfile.write(TEMPLATES.htmlSummary().format(
            # mod_file=os.path.abspath(globs['mod-file']),
            run_name=globs['run-name'],
            run_time=globs['startdatetimenice'],
            host_name=os.uname()[1],
            script_call=globs['call'],
            num_aln=str(globs['num-loci']),
            num_no_inf_loci=str(len(globs['no-inf-sites-loci'])),
            num_st_loci=str(globs['st-loci']),
            num_gt_loci=str(globs['gt-loci']),
            num_batches=str(globs['num-batches']),
            batch_size=str(globs['batch-size']),
            procs_per_batch=str(globs['procs-per-job']),
            num_st_batches=str(len(globs['st-batches'])),
            num_gt_batches=str(len(globs['gt-batches'])),
            num_jobs=str(globs['num-jobs']),
            num_spec=str(len(globs['tree-tips'])),
            num_targets=str(len(globs['targets'])),
            num_conserved=str(len(globs['conserved'])),
            num_outgroups=str(len(globs['outgroup'])),
            log_file=globs['logfilename'],
            aln_stats_file=globs['alnstatsfile'],
            snakemake_cmd=globs['smk-cmd'],
            theta_comment_start=theta_comment_start,
            theta_comment_end=theta_comment_end,
            coal_tree_comment_start=coal_tree_comment_start,
            coal_tree_file=globs['coal-tree-file'],
            coal_tree_comment_end=coal_tree_comment_end,
            input_tree_plot=os.path.join("plots", globs['input-tree-plot-file']),
            avg_aln_len=str(round(globs['avg-aln-len'], 3)),
            median_aln_len=str(round(globs['med-aln-len'], 3)),
            avg_seq_len_nogap=str(round(globs['avg-nogap-seq-len'], 3)),
            med_seq_len_nogap=str(round(globs['med-nogap-seq-len'], 3)),
            aln_len_hist=os.path.join("plots", globs['aln-len-plot-file']),
            seq_len_hist=os.path.join("plots", globs['seq-len-plot-file']),
            informative_sites_hist=os.path.join("plots", globs['inf-sites-plot-file']),
            informative_sites_frac_hist=os.path.join("plots", globs['inf-sites-frac-plot-file']),
            variable_informative_sites_plot=os.path.join("plots", globs['var-inf-sites-plot-file']),
            comment_start=comment_start,
            comment_end=comment_end,
            avg_scf_hist=os.path.join("plots", globs['avg-scf-hist-file']),
            low_scf_hist=os.path.join("plots", globs['low-scf-hist-file']),
            scf_tree_plot=os.path.join("plots", globs['scf-tree-plot-file']),
            bl_scf_plot=os.path.join("plots", globs['bl-scf-plot-file']),
            date_time=PC.getFooterDateTime(),
        ))
    # Write the HTML summary file using the template

    step_start_time = PC.report_step(globs, step, step_start_time, "Success");
    globs['html-summary-written'] = True;
    # Status update

    return globs;