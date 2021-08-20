############################################################
# For phyloacc site, 08.21
# This generates the file "people.html"
############################################################

import sys, os
sys.path.append('..')
import lib.read_chunks as RC

######################
# HTML template
######################

html_template = """
<!doctype html>
    {head}

<body>
	<div class="row" id="top_grid">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_header">PhyloAcc</div>
		<div class="col-2-24" id="margin"></div>
	</div>

	{nav}

	<div class="sep_div"></div>

	<div class="row">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_col">
            <h2>The team:</h2>
			<ul>
				<li><a href="https://edwards.oeb.harvard.edu/people/scott-v-edwards" target="_blank">Scott Edwards</a></li>
				<li><a href="http://www.stat.tsinghua.edu.cn/en/teambuilder/faculty/zhirui-hu/" target="_blank">Zhirhui Hu</li>
				<li><a href="https://statistics.fas.harvard.edu/people/taehee-lee" target="_blank">Taehee Lee</a></li>
				<li><a href="http://sites.fas.harvard.edu/~junliu/" target="_blank">Jun Liu</a></li>
				<li><a href="https://scholar.harvard.edu/tsackton/home" target="_blank">Tim Sackton</a></li>
				<li><a href="https://gwct.github.io/" target="_blank">Gregg Thomas</a></li>
				<li><a href="https://statistics.fas.harvard.edu/people/han-yan" target="_blank">Han Yan</a></li>
			</ul>
		</div>
        <div class="col-2-24" id="margin"></div>
	</div>

	<div class="sep_div"></div>

    {footer}
</body>
"""

######################
# Main block
######################
pagefile = "people.html";
print("Generating " + pagefile + "...");
title = "PhyloAcc - people"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));