# Python script to add PhysiPKPD to a PhysiCell project

import sys
import argparse
import requests
import os
from zipfile import ZipFile 

from helpers_for_get_and_add import * # contains append_suffix, common_flags, extract_from_url, etc

parser = argparse.ArgumentParser()

parser.add_argument('DIR', type=str, help='PhysiCell project directory to add PhysiPKPD capabilities')
common_flags(parser)

args = parser.parse_args()
DIR = args.DIR

if os.path.exists(DIR) is False:
    print(f"The target project directory {DIR} cannot be found...\n\tCheck for typos and try again.")
    exit(-1)

# Get PhysiPKPD stuff
print("----------------------")
print("Now getting PhysiPKPD...")

get_pkpd(args, DIR)
update_physicell_files(DIR, project_loaded=True)

# Get studio stuff
dir_name = os.path.abspath(DIR).split(os.sep)[-1] # protects aginst using the current working directory "."
studio_dir = append_suffix(f"{DIR}/../{dir_name}-studio")
USE_STUDIO = args.studio
get_studio_pkpd(USE_STUDIO, studio_dir)
if USE_STUDIO:
    makefile = f"{DIR}/Makefile"
    with open(makefile, "r") as f:
        for line in f:
            if "PROGRAM_NAME :=" in line:
                project_name = line.split(" ")[-1].strip('\n')

# Print advice
i=1
print(f"You are all set!")
print(f"\t{i}. Move into your new project folder:")
i+=1
print(f"\t\tcd {DIR}")
print(f"\t2. (Re-)compile the project:")
i+=1
print(f"\t\tmake -j 8")
if args.studio:
    print(f"\t{i}. Edit them with studio. Make sure your config file is correct:")
    i+=1
    print(f"\t     For a sample project:")
    print(f"\t\t(MacOS/Unix)\tpython {studio_dir}/bin/studio.py -c ./config/PhysiCell_settings.xml -e {project_name} --pkpd")
    print(f"\t\t(Windows)\tpython {studio_dir}\\bin\studio.py -c .\config\PhysiCell_settings.xml -e {project_name}.exe --pkpd")
    print(f"\t     For a template project:")
    print(f"\t\t(MacOS/Unix)\tpython {studio_dir}/bin/studio.py -c ./config/PhysiCell_settings.xml -e {project_name} --pkpd")
    print(f"\t\t(Windows)\tpython {studio_dir}\\bin\studio.py -c .\config\PhysiCell_settings.xml -e {project_name}.exe --pkpd")
else:
    print(f"\t{i}. Edit according to https://github.com/drbergman/PhysiPKPD/releases/latest")
    i+=1
    print(f"Consider passing in the --studio flag next time to get studio set up for easier editing")
print_generic_pkpd_advice(i)
