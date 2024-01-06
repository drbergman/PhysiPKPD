# Python script to set up a PhysiCell project using PhysiPKPD

import sys
import argparse
import requests
import os
from zipfile import ZipFile 

from append_suffix import append_suffix
from append_suffix import common_flags
from append_suffix import extract_from_url

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

PKPD_TAG = args.pkpd_tag
PKPD_BRANCH = args.pkpd_branch

PKPD_USE_BRANCH = PKPD_BRANCH is not None
USE_LATEST = (PKPD_USE_BRANCH is False) and (PKPD_TAG is None)
USE_TAG = (USE_LATEST is False) and (PKPD_USE_BRANCH is False)

remote_url = None
if USE_LATEST == True:
    response = requests.get(f"https://api.github.com/repos/drbergman/PhysiPKPD/releases/latest")
    remote_url =  response.json()["zipball_url"]
elif USE_TAG:
    response = requests.get(f"https://api.github.com/repos/drbergman/PhysiPKPD/tags")
    if PKPD_TAG.startswith('v') is False:
        PKPD_TAG = f"v{PKPD_TAG}"
    for tag in response.json():
        if tag["name"] == PKPD_TAG:
            remote_url = tag["zipball_url"]
            break
    local_file = f"PhysiPKPD-{PKPD_TAG}"
    # remote_url = 'https://github.com/drbergman/PhysiPKPD/releases/download/' + PKPD_TAG + '/' + PKPD_TAG + '.zip'
else:
    local_file = f"PhysiPKPD-{PKPD_BRANCH}"
    remote_url = 'https://github.com/drbergman/PhysiPKPD/archive/refs/heads/' + PKPD_BRANCH + '.zip'

print(f"PKPD remote_url = {remote_url}")
target_dir = append_suffix("PhysiPKPD-TEMP")
extract_from_url(remote_url, target_dir)
# data = requests.get(remote_url)
# local_file = append_suffix(local_file,'.zip')
# with open(local_file, 'wb')as file:
#    file.write(data.content)

# temp_dir = append_suffix(target_dir)
# with ZipFile(local_file, 'r') as zObject: 
#     zObject.extractall(path=temp_dir)

print(f"extracted PhysiPKPD to {target_dir}")

import shutil

os.rename(f"{target_dir}/addons/PhysiPKPD",f"{DIR}/addons/PhysiPKPD")
os.removedirs(f"{target_dir}/addons")
os.rename(f"{target_dir}/sample_projects_physipkpd",f"{DIR}/sample_projects_physipkpd")
os.rename(f"{target_dir}/LICENSE",f"{DIR}/addons/PhysiPKPD/LICENSE")
os.rename(f"{target_dir}/README.md",f"{DIR}/addons/PhysiPKPD/README.md")
shutil.rmtree(target_dir)
print(f"Moved PhysiPKPD files to {DIR}. Deleted {target_dir}")

# Update Makefile
print("----------------------")
print(f"Now updating {DIR}/sample_projects/Makefile-default and {DIR}/Makefile to be ready to make PhysiPKPD samples and projects.")

source_file = open(f'{DIR}/addons/PhysiPKPD/Makefile-PhysiPKPD-Samples.txt', "r")
with open(f'{DIR}/sample_projects/Makefile-default', 'a') as f:
    f.write("\n")
    shutil.copyfileobj(source_file, f)

os.rename(f'{DIR}/Makefile',f'{DIR}/Makefile-backup')
shutil.copyfile(f'{DIR}/sample_projects/Makefile-default', f'{DIR}/Makefile')

print(f"Updated Makefile to be ready for PhysiPKPD samples")

# Get studio stuff
studio_dir = None
USE_STUDIO = args.studio
if USE_STUDIO:
    print("----------------------")
    print("Now getting studio with physipkpd...")
    response = requests.get(f"https://api.github.com/repos/drbergman/PhysiCell-Studio/releases")
    max_pkpd_version = [0,0,0]
    remote_url = None
    def ver_comp(old,new):
        if new[0] > old[0]:
            return True
        elif new[0] < old[0]:
            return False
        elif len(old) == 1 and len(new) == 1:
            print("Two idential versions of studio found??")
            exit()
        else:
            return ver_comp(old[1:],new[1:])
    for release in response.json():
        if "pkpd" in release["tag_name"]:
            print(f"tag_name = {release['tag_name']}")
            version_ind = release["tag_name"].find('-v') + 2
            version_str = release["tag_name"][version_ind:]
            print(f"version_str = {version_str}")
            version = [int(x) for x in version_str.split('.')]
            if ver_comp(max_pkpd_version, version):
                max_pkpd_version = version
                remote_url = release["zipball_url"]
    if remote_url is None:
        print("No studio-pkpd release found???")
        USE_STUDIO = False
    else:
        print(f"studio-pkpd remote_url = {remote_url}")
        studio_dir = append_suffix(f"{DIR}-studio")
        extract_from_url(remote_url, studio_dir)

        # data = requests.get(remote_url)
        # local_file = append_suffix("studio-pkpd-TEMP",'.zip')
        # with open(local_file, 'wb')as file:
        #     file.write(data.content)

        # studio_dir_temp = append_suffix(f"{DIR}-studio-TEMP")
        # with ZipFile(local_file, 'r') as zObject: 
        #     zObject.extractall(path=studio_dir_temp)
        # folder_name = os.listdir(f"./{studio_dir_temp}")[0]
        # os.rename(f"{studio_dir_temp}/" + folder_name, studio_dir)

        # os.removedirs(studio_dir_temp)

print(f"You are all set!")
print(f"\t1. Move into your new project folder:")
print(f"\t\tcd {DIR}")
print(f"\t2. Make a sample project or a template project and make your PhysiPKPD model! Examples:")
print(f"\t\tmake pkpd-proliferation-sample\n\t\tmake pkpd-apoptosis-sample\n\t\tmake pkpd-template\n\n")
print(f"\t3. Compile the project:")
print(f"\t\tmake -j 8")
if USE_STUDIO:
    print(f"\t4. Edit them with studio:")
    print(f"\t  For a sample project:")
    print(f"\t\t(MacOS/Unix) python ../{studio_dir}/bin/studio.py -c ./config/pkpd_model.xml -e pkpd_sample --pkpd")
    print(f"\t\t(Windows) python ..\{studio_dir}\bin\studio.py -c .\config\pkpd_model.xml -e pkpd_sample.exe --pkpd\n")
    print(f"\t  For a template project:")
    print(f"\t\t(MacOS/Unix) python ../{studio_dir}/bin/studio.py -c ./config/pkpd_model.xml -e pkpd_project --pkpd")
    print(f"\t\t(Windows) python ..\{studio_dir}\bin\studio.py -c .\config\pkpd_model.xml -e pkpd_project.exe --pkpd\n\n")
else:
    print(f"\t4a. Run the samples with")
    print("\t\t(MaxOS/Unix) ./pkpd_sample ./config/pkpd_model.xml")
    print(f"\t\t(Windows) pkpd_sample.exe .\config\pkpd_model.xml")
    print(f"\t4b. Run the template projects with")
    print(f"\t\t (MaxOS/Unix) ./pkpd_project ./config/pkpd_model.xml")
    print(f"\t\t (Windows) pkpd_project.exe .\config\pkpd_model.xml")
print(f"\t5. Make sure to do the following:")
print(f"\t\t* Activate Dirichlet conditions for any PK substrate (that's how they enter the microenvironnent)")
print(f"\t\t* Set nonzero uptake rates for every (cell type, substrate) pairing that has a PD model (that's how damage is accumulated)")
print(f"\t\t* Add rules for how `custom:S_damage` (S = substrate name) affects target cell types; enable rules")
print(f"\t6. Consider using `damage_coloring` to assess your PD model parameters as you build the model.")
print(f"\t     In custom_modules/custom.cpp, replace my_coloring_function with\n")
print("std::vector<std::string> my_coloring_function(Cell *pC)\n{ return damage_coloring(pC); }\n\n")