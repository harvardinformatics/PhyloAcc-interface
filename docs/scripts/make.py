import sys, os, argparse

print()
print("###### Build site pages ######");
print("PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])))
print("# Script call: " + " ".join(sys.argv) + "\n----------");

parser = argparse.ArgumentParser(description="Gets stats from a bunch of abyss assemblies.");
parser.add_argument("--all", dest="all", help="Build all pages", action="store_true", default=False);
parser.add_argument("--index", dest="index", help="Without --all: build index.html. With --all: exlude index.html", action="store_true", default=False);
parser.add_argument("--cf", dest="cf", help="Without --all: build cf.html. With --all: exlude cf.html", action="store_true", default=False);
parser.add_argument("--mammals", dest="mammals", help="Without --all: build mammals.html. With --all: exlude mammals.html", action="store_true", default=False);
parser.add_argument("--sims", dest="sims", help="Without --all: build sim_trees.html. With --all: exlude sim_trees.html", action="store_true", default=False);
parser.add_argument("--people", dest="people", help="Without --all: build people.html. With --all: exlude people.html", action="store_true", default=False);
parser.add_argument("--links", dest="links", help="Without --all: build links.html. With --all: exlude links.html", action="store_true", default=False);
args = parser.parse_args();
# Input options.

#cwd = os.getcwd();
os.chdir("generators");

pages = {
    'index' : args.index,
    'cf' : args.cf,
    'mammals' : args.mammals,
    'sims' : args.sims,
    'people' : args.people,
    'links' : args.links
}

if args.all:
    pages = { page : False if pages[page] == True else True for page in pages };

if pages['index']:
    os.system("python index_generator.py");

if pages['cf']:
    os.system("Rscript cf_generator.r");

if pages['mammals']:
    os.system("Rscript mammals_generator.r");

if pages['sims']:
    os.system("Rscript sims_generator.r");

if pages['people']:
    os.system("python people_generator.py");

if pages['links']:
    os.system("python links_generator.py");
    
print("----------\nDone!");


