#############################################################################
# Functions to read sequences for PhyloAcc
# Gregg Thomas
#############################################################################

import sys
import os
import gzip
import lib.core as PC
import multiprocessing as mp
from itertools import groupby

############################################################################# 

def readFasta(filename, globs):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 

    if globs['seq-compression'] == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif globs['seq-compression'] == "none":
        file_stream = open(filename); 
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line[0] == ">"));
        readstr = lambda s : s.strip();
    # Read the lines of the file depending on the compression level
    # file_stream opens the file as an iterable
    # groupby takes an iterable (file_stream) and a function that indicates the key of the group. It iterates over
    # the iterable and when it encounters a key, it groups all following items to it (perfect for FASTA files).
    # fa_iter is a generator object holding all the grouped iterators.
    # readstr is a function that changes depending on compression level -- for compressed files we also need to decode
    # each string in the iterators below.

    seqdict = {};
    # A dictionary of sequences:
    # <sequence id/header> : <sequence>

    for header_obj in fa_iter:
        #print(header_obj)
        header = readstr(header_obj.__next__());
        # The header object is an iterator. This gets the string.

        curkey = header[1:];
        # This removes the ">" character from the header string to act as the key in seqdict

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        if curkey in globs['tree-tips']:
            seqdict[curkey] = seq;
        # Save the sequence in seqdict if it belongs to a tip branch in the input tree  

    return seqdict;

#############################################################################

def readBed(filename, globs):
# A function to read a bed file and store relevant info in a dict

    if globs['bed-compression'] == "gz":
        reader = gzip.open;
        readline = lambda l : l.decode().strip().split("\t");
    elif globs['bed-compression'] == "none":
        reader = open;
        readline = lambda l : l.strip().split("\t");
    # Read the lines of the file depending on the compression level

    bed_coords = {};
    # The bed information stored as a dict:
    # <locus id> : { <scaffold id>, <start coord>, <end coord> }

    first = True;
    # Flag to indicate the first line

    bed_id_field = True;
    # Flag to indicate whether the bed file contains locus IDs in the fourth column
    # Assume True, but check the first line for it below

    for line in reader(filename):
        line = readline(line);
        # Parse the current line into a list

        if first:
            if len(line) < 4:
                locus_id = 1;
                bed_id_field = False;
            # If there is no ID column with locus IDs, just set a counter as the ID here and
            # Set the id field flag to False
            first = False;
        # Check for an ID column in the first line

        if bed_id_field:
            cur_locus = line[3];
            if globs['id-file'] and cur_locus not in globs['locus-ids']:
                continue;
            # If the user provided and ID file to run a subset of loci in the input, check for
            # the current ID in that list here. If it's not in that list, do not save the locus
            # into the dictionary.
        # Get the ID from the current line
        else:
            cur_locus = str(locus_id);
            locus_id += 1;
        # If there is no ID column, set the ID here and increment

        scaffold, start, end = line[0], line[1], line[2];
        bed_coords[cur_locus] = { "scaff" : scaffold, "start" : int(start), "end" : int(end) };
        # Get info from the current line and save it in the dictionary

    return bed_coords;

#############################################################################

def partitionSeqs(concat_seqs, bed_coords):
# A function that takes a concatenated input sequence and coordinates from a bed file
# stored in a dictionary and separates the alignments by those coordinates

    alns = {};
    # The dictionary of individual alignments to return:
    # <locus id> : { <sequence id> : <sequence> }

    for locus in bed_coords:
    # Get every locus in the input coordinates
        cur_aln = {};
        # The dictionary to store the current locus alignment:
        # <sequence id> : <sequence>

        cur_start = bed_coords[locus]['start'];
        cur_end = bed_coords[locus]['end'];
        # Get the coordinates for the current locus

        for header in concat_seqs:
            cur_aln[header] = concat_seqs[header][cur_start:cur_end];
        # Get every sequence at the coordinates for the curent locus

        alns[locus] = cur_aln;
        # Save the current alignment to the alns dict

    return alns;

#############################################################################

def readSeq(globs):

    if globs['aln-file']:
        step = "Detecting compression of seq file";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        globs['seq-compression'] = PC.detectCompression(globs['aln-file']);
        if globs['seq-compression'] == "none":
            step_start_time = PC.report_step(globs, step, step_start_time, "Success: No compression detected");
        else:
            step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + globs['seq-compression'] + " detected");
        # Detect the compression of the input sequence file

        step = "Reading input FASTA";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        globs['in-seqs'] = readFasta(globs['aln-file'], globs);
        step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(globs['in-seqs'])) + " seqs read");
        # Read the input sequence file

        if globs['id-file']:
            step = "Reading locus IDs";
            step_start_time = PC.report_step(globs, step, False, "In progress...");
            globs['locus-ids'] = open(globs['id-file'], "r").read().split("\n");
            globs['locus-ids'] = list(filter(None, globs['locus-ids']));
            # Removes empty strings from the list of IDs
            step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(globs['locus-ids'])) + " IDs read");
        # If an ID file is provided, read the locus IDs here

        step = "Detecting compression of bed file";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        globs['bed-compression'] = PC.detectCompression(globs['bed-file']);
        if globs['bed-compression'] == "none":
            step_start_time = PC.report_step(globs, step, step_start_time, "Success: No compression detected");
        else:
            step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + globs['bed-compression'] + " detected");
        # Detect the compression of the input bed file

        step = "Reading input bed file";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        globs['in-bed'] = readBed(globs['bed-file'], globs);
        step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(globs['in-bed'])) + " loci read");
        # Read the input bed file with partition coordinates of the input alignment

        step = "Partitioning alignments by locus";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        globs['alns'] = partitionSeqs(globs['in-seqs'], globs['in-bed']);
        step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(globs['alns'])) + " alignments partitioned");
        # Separate the concatenated alignment to individual locus alignments based on the partitions in the bed file

        ##
        # step = "Writing partitioned sequences";
        # step_start_time = PC.report_step(globs, step, False, "In progress...");
        # written = 0;
        # for aln in globs['alns']:
        #     outdir = "simu_500_200_diffr_2-1";
        #     outfile = os.path.join(outdir, aln + ".fa");
        #     with open(outfile, "w") as of:
        #         for header in globs['alns'][aln]:
        #             of.write(">" + header + "\n");
        #             of.write(globs['alns'][aln][header] + "\n");
        #     written += 1;
        # step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(written) + " alignments written");
        # Chunk of code to write out the sequences in a concatenated file to individual files by locus -- for development
        ###

    # Read sequences if input is concatenated alignment + partitions in bed file
    #######################

    elif globs['aln-dir']:
        
        step = "Getting FASTA files in input dir";
        step_start_time = PC.report_step(globs, step, False, "In progress...");        
        aln_files = [ os.path.join(globs['aln-dir'], f) for f in os.listdir(globs['aln-dir']) if f.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")) ];
        step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(aln_files)) + " FASTA files found");
        
        step = "Detecting compression of seq files";
        step_start_time = PC.report_step(globs, step, False, "In progress...");
        globs['input-compression'] = PC.detectCompression(aln_files[0]);
        if globs['input-compression'] == "none":
            step_start_time = PC.report_step(globs, step, step_start_time, "Success: No compression detected");
        else:
            step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + globs['input-compression'] + " detected");
        # Detect the compression of the input sequence file

        step = "Reading input FASTA files";
        step_start_time = PC.report_step(globs, step, False, "In progress...");

        globs['alns'] = {};
        # The dictionary of individual alignments to return:
        # <locus id> : { <sequence id> : <sequence> }

        for f in aln_files:
            locus_id = os.path.splitext(os.path.basename(f))[0];
            # Get the locus ID from the file name

            cur_aln = readFasta(f, globs);
            # Read the current file as a FASTA file
            
            globs['alns'][locus_id] = cur_aln;
            # Add the current alignment to the main aln dict

        step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(globs['alns'])) + " files read");

    # Read sequences if input is a directory of alignment files
    #######################

    return globs;

#############################################################################

def locusAlnStats(locus_item):
    locus, aln, skip_chars = locus_item;
    # Unpack the data for the current locus

    num_seqs = len(aln);
    aln_len = len(aln[list(aln.keys())[0]]);

    cur_stats = { 'num-seqs' : num_seqs, 'length' : aln_len, 'variable-sites' : 0, 
                                        'informative-sites' : 0, 'num-sites-w-gap' : 0, 'num-sites-half-gap' : 0,
                                        'num-seqs-half-gap' : 0 };
    # Initialize the stats dict for this locus

    half_aln_len = aln_len / 2;
    half_site_len = num_seqs / 2;
    # Compute half the alignment length and half the site length for the current locus.
    # Used to see if seqs and sites are half-gap or more

    for j in range(aln_len):
        site = "";
        for seq in aln:
            if aln[seq][j].count("-") >= half_aln_len:
                cur_stats['num-seqs-half-gap'] += 1;
            # Count whether the current sequence is more than half gap

            site += aln[seq][j];
        # Get each allele from each sequence as the site

        allele_counts = { allele : site.count(allele) for allele in site if allele not in skip_chars };
        # Count the occurrence of each allele in the site

        if len(allele_counts) > 1:
            cur_stats['variable-sites'] += 1;
            # If there is more than one allele in the site, it is variable

            multi_allele_counts = [ allele for allele in allele_counts if allele_counts[allele] >= 2 ];
            # Count the number of allele present in at least 2 species

            if len(multi_allele_counts) >= 2:
                cur_stats['informative-sites'] += 1;
            # If 2 or more alleles are present in 2 or more species, this site is informative

        if "-" in site:
            cur_stats['num-sites-w-gap'] += 1;

            if site.count("-") >= half_site_len:
                cur_stats['num-sites-half-gap'] += 1;
        # Count whether this site contains a gap and, if so, whether more than half the sequences are a gap
    ## End site loop

    return locus, cur_stats;

#############################################################################

def alnStats(globs):
    step = "Calculating alignment stats";
    step_start_time = PC.report_step(globs, step, False, "In progress...");

    with globs['aln-pool'] as pool:
        for result in pool.imap(locusAlnStats, ((locus, globs['alns'][locus], globs['skip-chars']) for locus in globs['alns'])):
        # Loop over every locus in parallel to calculate stats
        # Have to do it this way so it doesn't terminate the pool for sCF calculations

            globs['aln-stats'][result[0]] = result[1];
            # Unpack the current result

    step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(globs['aln-stats'])) + " alignments processed");

    return globs;

#############################################################################

# read_seq = True;
#     # We only want to read sequences that are tips in the input tree. This flag
#     # keeps track of that

#     # for line in reader(filename):
#     #     line = readline(line);

#     tot_lines = 22336350
#     tot_lines_str = "22336350";
#     i = 0;

#     #for line in seqlines:
#     with open(filename) as seqfile:
#         for line in seqfile:
#             if i % 1000 == 0:
#                 print( str(i) + " / " + tot_lines_str + " (" + str((i / tot_lines) * 100) + ")")
#             i += 1;
#             line = line.strip();
#             #line = readline(line);

#             # if line in ["\n", ""]:
#             #     continue;
#             # Skip empty lines in the input

#             if line[0] == ">":
#                 curkey = line.replace(">", "");
#                 print();
#                 print(curkey);
#                 seqdict[curkey] = "";
#                 continue;

#                 # print()
#                 # print(curkey);
#                 # if curkey in globs['tree-tips']:
#                 #     read_seq = True;
#                 #     seqdict[curkey] = "";
#                 # else:
#                 #     read_seq = False;
#             # If the current line is a header, set that header as the current key and check
#             # that it is a tip in the input tree
#             #elif read_seq:
#             seqdict[curkey] += line;
#             # Otherwise, if the current header is a tip in the input tree the line contains 
#             # sequence that should be added to the current header's sequence






# for locus in globs['alns']:
#         globs['aln-stats'][locus] = { 'num-seqs' : len(globs['alns'][locus]), 'length' : "NA", 'variable-sites' : 0, 
#                                         'informative-sites' : 0, 'num-sites-w-gap' : 0, 'num-sites-half-gap' : 0};
#         # Calculate some basic alignment statistics for each locus
        
#         globs['aln-stats'][locus]['length'] = len(globs['alns'][locus][list(globs['alns'][locus].keys())[0]]);
#         # Get the length of the alignment from the first sequence

#         for j in range(globs['aln-stats'][locus]['length']):
#             site = "";
#             for seq in globs['alns'][locus]:
#                 site += globs['alns'][locus][seq][j];
#             # Get each allele from each sequence as the site

#             allele_counts = { allele : site.count(allele) for allele in site if allele not in globs['skip-chars'] };
#             # Count the occurrence of each allele in the site

#             if len(allele_counts) > 1:
#                 globs['aln-stats'][locus]['variable-sites'] += 1;
#                 # If there is more than one allele in the site, it is variable

#                 multi_allele_counts = [ allele for allele in allele_counts if allele_counts[allele] >= 2 ];
#                 # Count the number of allele present in at least 2 species

#                 if len(multi_allele_counts) >= 2:
#                     globs['aln-stats'][locus]['informative-sites'] += 1;
#                 # If 2 or more alleles are present in 2 or more species, this site is informative

#             if "-" in site:
#                 globs['aln-stats'][locus]['num-sites-w-gap'] += 1;

#                 if site.count("-") >= len(site) / 2:
#                     globs['aln-stats'][locus]['num-sites-half-gap'] += 1;
#             # Count whether this site contains a gap and, if so, whether more than half the sequences are a gap
#         ## End site loop
#     ## End locus loop