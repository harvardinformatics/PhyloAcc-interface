#!/usr/bin/env python3
#############################################################################
# A script to combine and summarize the outputs of a batched PhyloAcc run 
# with the interface
#
# Gregg Thomas
# Fall 2021
#############################################################################

import sys
import os
import lib.core as CORE
import phyloacc_lib.post_params as params
import phyloacc_lib.post_opt_parse as OP

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.

    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# PhyloAcc interface version " + globs['interface-version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual PhyloAcc version for this, and not just the interface version.

    print("#");
    print("# " + "=" * 150);

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    ####################

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    step = "Preparing output files";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    file_suffixes = ["_elem_lik.txt", "_M0_elem_Z.txt", "_M1_elem_Z.txt", "_M2_elem_Z.txt", "_rate_postZ_M0.txt", "_rate_postZ_M1.txt", "_rate_postZ_M2.txt"];
    # The files created by the various PhyloAcc runs that we care about have these suffixes
    # _elem_lik.txt MUST be first to get the locus IDs and associated numbers for each batch
        
    outfiles = [];
    # Combined output files will be created and stored in this list

    for f in file_suffixes:
        cur_outfile = os.path.join(globs['outdir'], f[1:]);
        open(cur_outfile, "w").close()
        outfiles.append(cur_outfile);
    # Create the combined output file for each suffix and add it to the list

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: files created");
    # Status update

    ####################

    step = "Combining batch outputs";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    batch_dirs = os.listdir(globs['phyloacc-out-dir']);
    # The directories within the PhyloAcc job directory, 1 for each batch

    id_keys = {};
    # Most output files only have locus numbers assigned by PhyloAcc. This dict stores the association with the actual
    # locus ID found in the _elem_lik.txt file
    # key:value format: <batch> : { <locus number> : <locus id> }

    first_file = True;
    # The first file flag - must be _elem_lik.txt

    for f in range(len(file_suffixes)):
    ## Loop through all file types. Not sure if I mentioned this before, but _elem_lik.txt MUST be first

        first_batch = True;
        # A first batch flag - we only write the headers in the output file during the first batch

        with open(outfiles[f], "a") as outfile:
        ## Open the current output file

            for batch_dir in batch_dirs:
            ## Go through every batch for the current file

                batch = batch_dir[0:batch_dir.index("-")];
                # Parse the batch string (just a number)

                if first_file:
                    id_keys[batch] = {};
                # If this is the first file (_elem_lik.txt), initialize the sub-dict in id_keys

                full_batch_dir = os.path.join(globs['phyloacc-out-dir'], batch_dir);
                cur_batch_file = os.path.join(full_batch_dir, batch + file_suffixes[f]);
                # Get the full batch directory and current batch file

                first_line = True;
                # The first line flag - the first line of each file contains the headers

                for line in open(cur_batch_file):
                ## Loop through every line in the current batch file

                    line = line.strip().split("\t");
                    # Parse the line into a list

                    if first_line:
                    ## If this is the first line (headers), we need to do some stuff, or just skip

                        if first_batch:
                        ## If this is also the first batch, we need to write the headers to the output file

                            cur_headers = line;
                            if first_file:
                                cur_headers = cur_headers[1:];
                            # Get the headers, and if this is the first file (_elem_lik.txt), remove the "No." header

                            cur_headers[0] = "Locus ID"
                            # Replace the "No." header with the "Locus ID header"

                            outfile.write("\t".join(cur_headers) + "\n");
                            # Write out the headers
                        ## End first batch block

                        first_line = False;
                        continue;
                        # Skip to the next line
                    ## End first line block

                    if first_file:
                        id_keys[batch][line[0]] = line[1];
                        outline = line[1:];
                    # If this is the first file (_elem_lik.txt), save the locus number:locus id in the id_keys dict and
                    # remove the "No." entry from the line
                    else:
                        outline = [id_keys[batch][line[0]]] + line[1:];
                    # For all other files, replace the locus number with the locus ID based on the entry in id_keys

                    outfile.write("\t".join(outline) + "\n");
                    # Write the line to the output file
                ## End line loop

                first_batch = False;
                # Switch off the batch flag
            ## End batch loop

            first_file = False;
            # Switch off the file flat
        ## Close batch output file
    ## End file loop

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: files combined");
    # Status update

    ####################

    CORE.endProg(globs, interface=False);

#############################################################################