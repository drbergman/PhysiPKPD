# Python script to add PhysiPKPD to a PhysiCell project

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
parser.add_argument('DIR', type=str, help='PhysiCell project directory to add PhysiPKPD capabilities')

group = parser.add_mutually_exclusive_group()
group.add_argument('-pt','--pkpd_tag', default=None, help='the tag of the PhysiCell release to use')
group.add_argument('-pb','--pkpd_branch', default=None, help='use this branch of a PhysiCell fork rather than a release')

parser.add_argument('--studio',action='store_true', help='also downloads a copy of studio with physipkpd integration')

args = parser.parse_args()
DIR = args.DIR

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

os.rename(f"{temp_dir}/{folder_name}/addons/PhysiPKPD",f"{DIR}/addons/PhysiPKPD")
os.removedirs(f"{temp_dir}/{folder_name}/addons")
# os.rename(f"{temp_dir}/{folder_name}/sample_projects_physipkpd",f"{DIR}/sample_projects_physipkpd") # don't need to move sample projects in this case
os.rename(f"{temp_dir}/{folder_name}/LICENSE",f"{DIR}/addons/PhysiPKPD/LICENSE")
os.rename(f"{temp_dir}/{folder_name}/README.md",f"{DIR}/addons/PhysiPKPD/README.md")
shutil.rmtree(temp_dir)
print(f"Moved PhysiPKPD files to {DIR}. Deleted {temp_dir}")

source_file = f'{DIR}/addons/PhysiPKPD/Makefile-PhysiPKPD-Objects.txt'
with open(source_file, 'r') as sf:
    new_lines = sf.readlines()
print(new_lines)
exit()
with open(f'{DIR}/Makefile', 'r') as f:
    lines = f.readlines()

for line_number, line in enumerate(lines):
    if "PhysiCell_custom_module_OBJECTS :=" in line:
        break

with open(f'{DIR}/Makefile','w') as f:
    for i in range(line_number):
        print(lines[i], file=f)
    for new_line in new_lines:
        print(new_line, file=f)
    for i in range(i+1,len(lines)):
        print(lines[i], file=f)

with open(f'{DIR}/Makefile','r+') as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "PhysiCell_OBJECTS :=" in line:
        line = line.rstrip('\n') + ' $(PhysiPKPD_OBJECTS)\n'
        break
lines[i] = line

with open(f'{DIR}/Makefile','w') as f:
    f.writelines(lines)

print(f"Updated Makefile to be ready for PhysiPKPD")

def get_line_number(s, lines):
    for line_number, line in enumerate(lines):
        if s in line:
            return line_number
    print(f"{f.name} does not have a line containing {s}")
    
with open(f'{DIR}/main.cpp', 'r+') as f:
    lines = f.readlines()
    line_number = get_line_number("setup_tissue", lines)
    if line_number is not None:
        lines.insert(line_number+1, "\n\tsetup_pharmacodynamics();\n")  # new_string should end in a newline
        f.seek(0)  # readlines consumes the iterator, so we need to start over
        f.writelines(lines)  # No need to truncate as we are increasing filesize
        
with open(f'{DIR}/main.cpp', 'r+') as f:
    lines = f.readlines()
    line_number = get_line_number("microenvironment.simulate_diffusion_decay", lines)
    if line_number is not None:
        lines.insert(line_number+1, "\n\t\t\tPD_model( PhysiCell_globals.current_time );\n")  # new_string should end in a newline
        lines.insert(line_number, "\t\t\tPK_model( PhysiCell_globals.current_time );\n\n")  # new_string should end in a newline
        f.seek(0)  # readlines consumes the iterator, so we need to start over
        f.writelines(lines)  # No need to truncate as we are increasing filesize

with open(f"{DIR}/custom_modules/custom.h", "r+") as f:
    lines = f.readlines()
    line_number = get_line_number("#include", lines)
    lines.insert(line_number,'#include "../addons/PhysiPKPD/src/PhysiPKPD.h"\n')
    f.seek(0)
    f.writelines(lines)

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

    studio_dir = append_suffix(f"{DIR}-studio")
    studio_dir_temp = append_suffix(f"{DIR}-studio-TEMP")
    with ZipFile(local_file, 'r') as zObject: 
        zObject.extractall(path=studio_dir_temp)
    folder_name = os.listdir(f"./{studio_dir_temp}")[0]
    os.rename(f"{studio_dir_temp}/" + folder_name, studio_dir)

    os.removedirs(studio_dir_temp)

print(f"You are all set!")
print(f"\t1. Move into your project folder:")
print(f"\t\tcd {DIR}")
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
