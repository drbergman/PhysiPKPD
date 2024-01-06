# Python script to set up a PhysiCell project using PhysiPKPD

import sys
import argparse
import requests
import os
from zipfile import ZipFile 
import shutil

from helpers_for_get_and_add import * # contains append_suffix, common_flags, extract_from_url, etc

parser = argparse.ArgumentParser()

group = parser.add_mutually_exclusive_group()
group.add_argument('-t','--tag', default=None, help='the tag of the PhysiCell release to use')
group.add_argument('-b','--branch', default=None, help='use this branch of a PhysiCell fork rather than a release')

parser.add_argument('-o','--owner', default='MathCancer', help='the owner of the PhysiCell fork')
parser.add_argument('-d','--dir', default='PhysiPKPD_Project',help='target directory for new PhysiCell folder')

common_flags(parser)

args = parser.parse_args()

OWNER = args.owner
TAG = args.tag
BRANCH = args.branch
DIR = args.dir

USE_BRANCH = BRANCH is not None
USE_LATEST = (USE_BRANCH is False) and (TAG is None)
USE_TAG = (USE_LATEST is False) and (USE_BRANCH is False)

remote_url = None
if USE_LATEST:
    response = requests.get(f"https://api.github.com/repos/{OWNER}/PhysiCell/releases/latest")
    remote_url = response.json()["zipball_url"]
elif USE_TAG:
    response = requests.get(f"https://api.github.com/repos/{OWNER}/PhysiCell/tags")
    for tag in response.json():
        if tag["name"] == TAG:
            remote_url = tag["zipball_url"]
            break
else:
    response = requests.get(f"https://api.github.com/repos/{OWNER}/PhysiCell/branches")
    print(response.json())
    for branch in response.json():
        print(f"branch = {branch}")
        if branch["name"] == BRANCH:
            remote_url = f"https://api.github.com/repos/{OWNER}/PhysiCell/zipball/{BRANCH}"
            break

print(f"remote_url = {remote_url}")
DIR = append_suffix(DIR)
extract_from_url(remote_url, DIR)
print("unzipped to ",DIR)
    
# Get PhysiPKPD stuff
print("----------------------")
print("Now getting PhysiPKPD...")

get_pkpd(args, DIR)
update_physicell_files(DIR)

# Get studio stuff
studio_dir = append_suffix(f"{DIR}-studio")
USE_STUDIO = args.studio
get_studio_pkpd(USE_STUDIO, studio_dir)

# Print advice
i=1
print(f"You are all set!")
print(f"\t{i}. Move into your new project folder:")
i+=1
print(f"\t\tcd {DIR}")
print(f"\t{i}. Make a sample project or a template project and make your PhysiPKPD model! Examples:")
i+=1
print(f"\t\tmake pkpd-proliferation-sample\n\t\tmake pkpd-apoptosis-sample\n\t\tmake pkpd-template")
print(f"\t{i}. Compile the project:")
i+=1
print(f"\t\tmake -j 8")
if USE_STUDIO:
    print(f"\t{i}. Edit them with studio:")
    i+=1
    print(f"\t  For a sample project:")
    print(f"\t\t(MacOS/Unix) python ../{studio_dir}/bin/studio.py -c ./config/pkpd_model.xml -e pkpd_sample --pkpd")
    print(f"\t\t(Windows) python ..\{studio_dir}\bin\studio.py -c .\config\pkpd_model.xml -e pkpd_sample.exe --pkpd\n")
    print(f"\t  For a template project:")
    print(f"\t\t(MacOS/Unix) python ../{studio_dir}/bin/studio.py -c ./config/pkpd_model.xml -e pkpd_project --pkpd")
    print(f"\t\t(Windows) python ..\{studio_dir}\bin\studio.py -c .\config\pkpd_model.xml -e pkpd_project.exe --pkpd\n\n")
else:
    print(f"\t{i}a. Run the samples with")
    print("\t\t(MaxOS/Unix) ./pkpd_sample ./config/pkpd_model.xml")
    print(f"\t\t(Windows) pkpd_sample.exe .\config\pkpd_model.xml")
    print(f"\t{i}b. Run the template projects with")
    i+=1
    print(f"\t\t (MaxOS/Unix) ./pkpd_project ./config/pkpd_model.xml")
    print(f"\t\t (Windows) pkpd_project.exe .\config\pkpd_model.xml")
print_generic_pkpd_advice(i)