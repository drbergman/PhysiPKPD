# Python script to set up a PhysiCell project using PhysiPKPD

import sys
import argparse
import requests
import os
from zipfile import ZipFile 

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

with ZipFile(local_file, 'r') as zObject: 
    zObject.extractall(path=DIR_NAME + "_temp")

suffix = ""
while True:
    try:
        os.rename(DIR_NAME + "_temp/" + zip_folder_name, DIR_NAME + str(suffix))
        break
    except:
        if suffix == "":
            suffix = 1
        else:
            suffix += 1
        print(f"Trying suffix {suffix} now...")

os.removedirs(DIR_NAME + "_temp")
DIR_NAME += str(suffix)
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
    print(f"response = {response}")
    release_name_str = response.json()["name"]
    print(release_name_str)
    print(release_name_str.split())
    vnum = release_name_str.split()[1]
    print("vnum =",vnum)  # e.g., vnum= 1.10.4
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
with open(local_file, 'wb')as file:
   file.write(data.content)


with ZipFile(local_file, 'r') as zObject: 
    zObject.extractall(path="PhysiPKPD-TEMP")

print(f"extracted PhysiPKPD to PhysiPKPD-TEMP")

folder_name = os.listdir("./PhysiPKPD-TEMP")[0]

import shutil

os.rename(f"PhysiPKPD-TEMP/{folder_name}/addons/PhysiPKPD",f"{DIR_NAME}/addons/PhysiPKPD")
os.removedirs(f"PhysiPKPD-TEMP/{folder_name}/addons")
os.rename(f"PhysiPKPD-TEMP/{folder_name}/sample_projects_physipkpd",f"{DIR_NAME}/sample_projects_physipkpd")
os.rename(f"PhysiPKPD-TEMP/{folder_name}/LICENSE",f"{DIR_NAME}/addons/PhysiPKPD/LICENSE")
os.rename(f"PhysiPKPD-TEMP/{folder_name}/README.md",f"{DIR_NAME}/addons/PhysiPKPD/README.md")
shutil.rmtree("PhysiPKPD-TEMP")
print(f"Moved PhysiPKPD files to {DIR_NAME}. Deleted PhysiPKPD-TEMP")

source_file = open(f'{DIR_NAME}/addons/PhysiPKPD/Makefile-PhysiPKPD_Addendum.txt', "r")
with open(f'{DIR_NAME}/sample_projects/Makefile-default', 'a') as f:
    f.write("\n")
    shutil.copyfileobj(source_file, f)


os.rename(f'{DIR_NAME}/Makefile',f'{DIR_NAME}/Makefile-backup')
shutil.copyfile(f'{DIR_NAME}/sample_projects/Makefile-default', f'{DIR_NAME}/Makefile')

print(f"Updated Makefile to be ready for PhysiPKPD samples")

print(f"You are all set!\n  Move into your new project folder:\n\tcd {DIR_NAME}\n  make a sample project or a template project and make your PhysiPKPD model!\n\n")
print(f"Examples:\n\tmake pkpd-proliferation-sample\n\tmake pkpd-apoptosis-sample\n\tmake pkpd-template\n\n")
print(f"Run the samples with\n\t./pkpd_sample ./config/pkpd_model.xml (MacOS/Linux)\n\tpkpd_sample.exe .\config\pkpd_model.xml (Windows)")
print(f"Run the template projects with\n\t./pkpd_project ./config/pkpd_model.xml (MacOS/Linux)\n\tpkpd_project.exe .\config\pkpd_model.xml (Windows)")
