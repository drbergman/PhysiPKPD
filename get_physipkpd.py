# Python script to set up a PhysiCell project using PhysiPKPD

import sys
import argparse
import requests
import os
from zipfile import ZipFile 

def append_suffix(f):
    suffix = ""
    while os.path.exists(f"{f}{suffix}"):
        if suffix == "":
            suffix = 1
        else:
            suffix += 1
    return f"{f}{suffix}"

parser = argparse.ArgumentParser()

group = parser.add_mutually_exclusive_group()
group.add_argument('-t','--tag', default=None, help='the tag of the PhysiCell release to use')
group.add_argument('-b','--branch', default=None, help='use this branch of a PhysiCell fork rather than a release')

parser.add_argument('-o','--owner', default='MathCancer', help='the owner of the PhysiCell fork')
parser.add_argument('-d','--dir', default='PhysiPKPD_Project',help='target directory for new PhysiCell folder')

group = parser.add_mutually_exclusive_group()
group.add_argument('-pt','--pkpd_tag', default=None, help='the tag of the PhysiCell release to use')
group.add_argument('-pb','--pkpd_branch', default=None, help='use this branch of a PhysiCell fork rather than a release')

args = parser.parse_args()

OWNER = args.owner
TAG = args.tag
BRANCH = args.branch
DIR_NAME = args.dir

USE_BRANCH = BRANCH is not None
USE_LATEST = (USE_BRANCH is False) and (TAG is None)
USE_TAG = (USE_LATEST is False) and (USE_BRANCH is False)

if USE_LATEST == True:
    response = requests.get(f"https://api.github.com/repos/{OWNER}/PhysiCell/releases/latest")
    release_name_str = response.json()["name"]
    print(release_name_str)
    print(release_name_str.split())
    vnum = release_name_str.split()[1]
    print("vnum =",vnum)  # e.g., vnum= 1.10.4
    zip_folder_name = 'PhysiCell'
    remote_url = f'https://github.com/{OWNER}/PhysiCell/releases/download/' + vnum + '/PhysiCell_V.' + vnum + '.zip'
elif USE_TAG:
    zip_folder_name = 'PhysiCell'
    remote_url = f'https://github.com/{OWNER}/PhysiCell/releases/download/' + TAG + '/PhysiCell_V.' + TAG + '.zip'
else:
    zip_folder_name = 'PhysiCell-' + BRANCH
    remote_url = f'https://github.com/{OWNER}/PhysiCell/archive/refs/heads/' + BRANCH + '.zip'

print("remote_url=",remote_url)
local_file = 'PhysiCell.zip'
data = requests.get(remote_url)
with open(local_file, 'wb')as file:
  file.write(data.content)

temp_dir = append_suffix(DIR_NAME + "-TEMP")
# suffix_temp_dir = 1
# while os.path.exists(temp_dir + str(suffix_temp_dir)):
#     suffix_temp_dir+=1
# temp_dir +=  str(suffix_temp_dir)

with ZipFile(local_file, 'r') as zObject: 
    zObject.extractall(path=temp_dir)

DIR_NAME = append_suffix(DIR_NAME)
os.rename(f"{temp_dir}/" + zip_folder_name, DIR_NAME)
# suffix = ""
# while True:
#     try:
#         os.rename(f"{temp_dir}/" + zip_folder_name, DIR_NAME + str(suffix))
#         break
#     except:
#         if suffix == "":
#             suffix = 1
#         else:
#             suffix += 1
#         print(f"Trying suffix {suffix} now...")

os.removedirs(temp_dir)
# DIR_NAME += str(suffix)
print("unzipped to ",DIR_NAME)
    
# Get PhysiPKPD stuff
print("----------------------")
print("Now getting PhysiPKPD...")

PKPD_TAG = args.pkpd_tag
PKPD_BRANCH = args.pkpd_branch

PKPD_USE_BRANCH = PKPD_BRANCH is not None
USE_LATEST = (PKPD_USE_BRANCH is False) and (PKPD_TAG is None)
USE_TAG = (USE_LATEST is False) and (PKPD_USE_BRANCH is False)

# physipkpd_branch = "development-rules-integration-studio-version"

if USE_LATEST == True:
    response = requests.get(f"https://api.github.com/repos/drbergman/PhysiPKPD/releases/latest")
    tag_name = response.json()["tag_name"]
    local_file = 'PhysiPKPD-LATEST.zip'
    remote_url =  response.json()["zipball_url"]
elif USE_TAG:
    local_file = f"PhysiPKPD-{PKPD_TAG}.zip"
    remote_url = 'https://github.com/drbergman/PhysiPKPD/releases/download/' + PKPD_TAG + '/' + PKPD_TAG + '.zip'
else:
    local_file = f"PhysiPKPD-{PKPD_BRANCH}.zip"
    remote_url = 'https://github.com/drbergman/PhysiPKPD/archive/refs/heads/' + PKPD_BRANCH + '.zip'

print("remote_url=",remote_url)
data = requests.get(remote_url)
local_file = append_suffix(local_file)
with open(local_file, 'wb')as file:
   file.write(data.content)

temp_dir = append_suffix("PhysiPKPD-TEMP")
with ZipFile(local_file, 'r') as zObject: 
    zObject.extractall(path=temp_dir)

print(f"extracted PhysiPKPD to {temp_dir}")

folder_name = os.listdir(f"./{temp_dir}")[0]

import shutil

os.rename(f"{temp_dir}/{folder_name}/addons/PhysiPKPD",f"{DIR_NAME}/addons/PhysiPKPD")
os.removedirs(f"{temp_dir}/{folder_name}/addons")
os.rename(f"{temp_dir}/{folder_name}/sample_projects_physipkpd",f"{DIR_NAME}/sample_projects_physipkpd")
os.rename(f"{temp_dir}/{folder_name}/LICENSE",f"{DIR_NAME}/addons/PhysiPKPD/LICENSE")
os.rename(f"{temp_dir}/{folder_name}/README.md",f"{DIR_NAME}/addons/PhysiPKPD/README.md")
shutil.rmtree(temp_dir)
print(f"Moved PhysiPKPD files to {DIR_NAME}. Deleted {temp_dir}")

source_file = open(f'{DIR_NAME}/addons/PhysiPKPD/Makefile-PhysiPKPD_Addendum.txt', "r")
with open(f'{DIR_NAME}/sample_projects/Makefile-default', 'a') as f:
    f.write("\n")
    shutil.copyfileobj(source_file, f)

os.rename(f'{DIR_NAME}/Makefile',f'{DIR_NAME}/Makefile-backup')
shutil.copyfile(f'{DIR_NAME}/sample_projects/Makefile-default', f'{DIR_NAME}/Makefile')

print(f"Updated Makefile to be ready for PhysiPKPD samples")

print(f"You are all set!\n\tMove into your new project folder:\n\t\tcd {DIR_NAME}\n\tmake a sample project or a template project and make your PhysiPKPD model!\n\n")
print(f"\tExamples:\n\t\tmake pkpd-proliferation-sample\n\t\tmake pkpd-apoptosis-sample\n\t\tmake pkpd-template\n\n")
print(f"\tRun the samples with\n\t\t./pkpd_sample ./config/pkpd_model.xml (MacOS/Linux)\n\t\tpkpd_sample.exe .\config\pkpd_model.xml (Windows)")
print(f"\tRun the template projects with\n\t\t./pkpd_project ./config/pkpd_model.xml (MacOS/Linux)\n\t\tpkpd_project.exe .\config\pkpd_model.xml (Windows)")
