############################################################
# For phyloacc site, 08.21
# This generates the file "index.html"
############################################################

import sys, os
sys.path.append(os.path.abspath('../lib/'))
import read_chunks as RC

######################
# HTML template
######################

html_template = """
<!doctype html>
    {head}

<body>
	<div class="row" id="top_grid">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_header">PhyloAcc interface</div>
		<div class="col-2-24" id="margin"></div>
	</div>

	{nav}

	<div class="sep_div"></div>

	<div class="row">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_col">
			<ul>
				<li><h3><a href="cf.html">Concordance factors to inform model selection</a></h3></li>
				<li><h3><a href="mammals.html">Accelerated mammal elements with gene tree model</a></h3></li>
			</ul>
		</div>
        <div class="col-2-24" id="margin"></div>
	</div>

	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>

    {footer}
</body>
"""

######################
# Main block
######################
pagefile = "index.html";
print("Generating " + pagefile + "...");
title = "PhyloAcc"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));