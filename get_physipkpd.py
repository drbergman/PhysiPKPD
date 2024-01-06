# Python script to set up a PhysiCell project using PhysiPKPD

import sys
import argparse
import requests
import os
from zipfile import ZipFile 

def append_suffix(f,ext=""):
    suffix = ""
    while os.path.exists(f"{f}{suffix}{ext}"):
        if suffix == "":
            suffix = 1
        else:
            suffix += 1
    return f"{f}{suffix}{ext}"

parser = argparse.ArgumentParser()

group = parser.add_mutually_exclusive_group()
group.add_argument('-t','--tag', default=None, help='the tag of the PhysiCell release to use')
group.add_argument('-b','--branch', default=None, help='use this branch of a PhysiCell fork rather than a release')

parser.add_argument('-o','--owner', default='MathCancer', help='the owner of the PhysiCell fork')
parser.add_argument('-d','--dir', default='PhysiPKPD_Project',help='target directory for new PhysiCell folder')

group = parser.add_mutually_exclusive_group()
group.add_argument('-pt','--pkpd_tag', default=None, help='the tag of the PhysiCell release to use')
group.add_argument('-pb','--pkpd_branch', default=None, help='use this branch of a PhysiCell fork rather than a release')

parser.add_argument('--studio',action='store_true', help='also downloads a copy of studio with physipkpd integration')

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

with ZipFile(local_file, 'r') as zObject: 
    zObject.extractall(path=temp_dir)

DIR_NAME = append_suffix(DIR_NAME)
os.rename(f"{temp_dir}/" + zip_folder_name, DIR_NAME)

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

if USE_LATEST == True:
    response = requests.get(f"https://api.github.com/repos/drbergman/PhysiPKPD/releases/latest")
    tag_name = response.json()["tag_name"]
    local_file = 'PhysiPKPD-LATEST'
    remote_url =  response.json()["zipball_url"]
elif USE_TAG:
    local_file = f"PhysiPKPD-{PKPD_TAG}"
    remote_url = 'https://github.com/drbergman/PhysiPKPD/releases/download/' + PKPD_TAG + '/' + PKPD_TAG + '.zip'
else:
    local_file = f"PhysiPKPD-{PKPD_BRANCH}"
    remote_url = 'https://github.com/drbergman/PhysiPKPD/archive/refs/heads/' + PKPD_BRANCH + '.zip'

print("remote_url=",remote_url)
data = requests.get(remote_url)
local_file = append_suffix(local_file,'.zip')
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

# Get studio stuff
studio_dir = None
if args.studio:
    print("----------------------")
    print("Now getting studio with physipkpd...")
    response = requests.get(f"https://api.github.com/repos/drbergman/PhysiCell-Studio/releases")
    max_pkpd_version = [0,0,0]
    zip_url = None
    def ver_comp(old,new):
        if new[0] > old[0]:
            return True
        elif new[0] < old[0]:
            return False
        elif len(old) == 1 and len(new) == 1:
            print("Two idential versions??")
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
                zip_url = release["zipball_url"]
    if zip_url is None:
        print("No studio-pkpd release found")
        exit()
    
    remote_url = zip_url
    print(f"remote_url = {remote_url}")
    data = requests.get(remote_url)
    local_file = append_suffix("studio-pkpd-TEMP",'.zip')
    with open(local_file, 'wb')as file:
        file.write(data.content)

    studio_dir = append_suffix(f"{DIR_NAME}-studio")
    studio_dir_temp = append_suffix(f"{DIR_NAME}-studio-TEMP")
    with ZipFile(local_file, 'r') as zObject: 
        zObject.extractall(path=studio_dir_temp)
    folder_name = os.listdir(f"./{studio_dir_temp}")[0]
    os.rename(f"{studio_dir_temp}/" + folder_name, studio_dir)

    os.removedirs(studio_dir_temp)

print(f"You are all set!")
print(f"\t1. Move into your new project folder:")
print(f"\t\tcd {DIR_NAME}")
print(f"\t2. Make a sample project or a template project and make your PhysiPKPD model! Examples:")
print(f"\t\tmake pkpd-proliferation-sample\n\t\tmake pkpd-apoptosis-sample\n\t\tmake pkpd-template\n\n")
print(f"\t3. Compile the project:")
print(f"\t\tmake -j 8")
if args.studio:
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