############################################################
# For murine web development, 11.17
# Functions to read static html chunks
############################################################

def readHead(title, active_url):
    page_name = active_url.replace(".html", "");
    headfile = "../html-chunks/head.html";
    return open(headfile, "r").read().replace("TMPTITLE", title).replace("PAGECSS", page_name);

def readNav(active_url):
    navfile = "../html-chunks/nav.html";
    navlines = open(navfile, "r").readlines()
    for x in range(len(navlines)):
        # if active_url != "index.html":
        if active_url in navlines[x]:
            navlines[x] = navlines[x].replace(active_url, "#");
            if 'id="nav_link_cell"' in navlines[x]:
                navlines[x] = navlines[x].replace('class="nav_link"', 'class="nav_link" id="active"');
                navlines[x] = navlines[x].replace('id="nav_link_cell"', 'id="nav_link_cell_active"');

        if active_url == "archive.html":
            if "dropdown" in navlines[x]:
                navlines[x] = navlines[x].replace('id="nav_link_cell"', 'id="nav_link_cell_active"');
            if "Updates" in navlines[x]:
                navlines[x] = navlines[x].replace('class="nav_link"', 'class="nav_link" id="active"');

        if "mobile_nav" in navlines[x]:
            break;
            
    return "".join(navlines);

def readLinks():
    linkfile = "../html-chunks/links.html";
    return open(linkfile, "r").read();

def readFooter():
    import time, subprocess
    from datetime import datetime
    footerfile = "../html-chunks/footer.html";
    now = datetime.now().strftime("%m/%d/%Y %H:%M:%S");
    zone = subprocess.check_output("date +%Z").decode().strip();
    return open(footerfile, "r").read().replace("DATETIME", now + " " + zone);