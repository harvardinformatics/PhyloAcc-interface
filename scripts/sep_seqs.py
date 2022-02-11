#############################################################################
# A script to separate the sequences for the mammal test data
# Gregg Thomas
#############################################################################

import os, sys
import seqparse as SEQ

#############################################################################

indir = "/n/holylfs05/LABS/informatics/Users/gthomas/PhyloAcc-interface-data/gt-test-mammals/chrome-data/";
outdir = "/n/holylfs05/LABS/informatics/Users/gthomas/PhyloAcc-interface-data/gt-test-mammals/seq-for-trees/";

chrome_dirs = [ d for d in os.listdir(indir) if d.startswith("chr") ];

for chrome in chrome_dirs:
    print(chrome);
    chrome_dir = os.path.join(indir, chrome);
    
    fa_file = os.path.join(chrome_dir, chrome + ".fasta");
    bed_file = os.path.join(chrome_dir, chrome + ".bed");

    if not os.path.isfile(fa_file) or not os.path.isfile(bed_file):
        print("SKIPPED");
        continue;

    cur_seqs = SEQ.fastaGetDict(fa_file);

    for line in open(bed_file):
        line = line.strip().split("\t");
        locus = line[0];
        start = int(line[1]);
        end = int(line[2]);

        #print(start, end);

        outfilename = os.path.join(outdir, chrome + "-" + locus + ".fa");
        #print(outfilename);
        with open(outfilename, "w") as outfile:
            for title in cur_seqs:
                cur_seq = cur_seqs[title][start:end];

                if set(cur_seq) in [{'*'}, {'-'}, {'N'}, {'*', '-'}, {'*', 'N'}, {'-', 'N'}, {'*', '-', 'N'}]:
                    continue;
                # For gene trees, skip sequences that are all gap

                outfile.write(title + "\n");
                outfile.write(cur_seq + "\n");
                #print(len(cur_seqs[title][start:end]));