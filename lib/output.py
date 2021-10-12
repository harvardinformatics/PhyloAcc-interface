#############################################################################
# Output functions for the PhyloAcc interface
# Gregg Thomas
#############################################################################

import os
import lib.core as PC

#############################################################################

def writeAlnStats(globs):

    step = "Writing: " + globs['alnstatsfile'];
    step_start_time = PC.report_step(globs, step, False, "In progress...");

    globs['alnstatsfile'] = os.path.join(globs['outdir'], globs['alnstatsfile']);
    loci_sorted = sorted(globs['aln-stats']);
    with open(globs['alnstatsfile'], "w") as outfile:
        first = True;
        for locus in loci_sorted:
            if first:
                keys = list(globs['aln-stats'][locus].keys())
                headers = ['locus'] + keys;
                outfile.write(",".join(headers) + "\n");
                first = False;
            outline = [locus] + [ str(globs['aln-stats'][locus][key]) for key in keys ];
            outfile.write(",".join(outline) + "\n");

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: align stats written");
    globs['aln-stats-written'] = True;

    return globs;

#############################################################################

def writeSCFStats(globs):

    step = "Writing: " + globs['scfstatsfile'];
    step_start_time = PC.report_step(globs, step, False, "In progress...");

    headers = ["node","variable-sites","decisive-sites","concordant-sites","quartet-scf-sum","num-quartets","total-quartets","avg-quartet-scf"];

    globs['scfstatsfile'] = os.path.join(globs['outdir'], globs['scfstatsfile']);
    with open(globs['scfstatsfile'], "w") as outfile:
        outfile.write(",".join(headers) + "\n");
        first = True;
        for node in globs['scf']:
            # if first:
            #     outfile.write(",".join(headers) + "\n");
            #     first = False;
            #     continue;

            globs['scf'][node]['num-quartets'] = len(globs['quartets'][node])

            outline = [node] + [ str(globs['scf'][node][header]) for header in headers if header != "node" ];
            outfile.write(",".join(outline) + "\n");

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: sCF stats written");
    globs['scf-stats-written'] = True;       


    step = "Writing: " + globs['scftreefile'];
    step_start_time = PC.report_step(globs, step, False, "In progress...");

    globs['scftreefile'] = os.path.join(globs['outdir'], globs['scftreefile']);
    labeled_scf_tree = globs['labeled-tree'];
    for node in globs['scf']:
        labeled_scf_tree = labeled_scf_tree.replace(node, node + "_" + str(round(globs['scf'][node]['avg-quartet-scf'], 2)));

    with open(globs['scftreefile'], "w") as outfile:
        outfile.write(labeled_scf_tree);
    step_start_time = PC.report_step(globs, step, step_start_time, "Success: sCF tree written");
    globs['scf-tree-written'] = True;   

    return globs;

## ALSO NEED TO WRITE TREE WITH CF, KEEPING LABELS IF PROVIDED

#############################################################################