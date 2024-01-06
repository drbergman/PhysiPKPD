import os
import argparse
import requests
from zipfile import ZipFile 

def append_suffix(f,ext=""):
    suffix = ""
    while os.path.exists(f"{f}{suffix}{ext}"):
        if suffix == "":
            suffix = 1
        else:
            suffix += 1
    return f"{f}{suffix}{ext}"

def common_flags(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-pt','--pkpd_tag', default=None, help='the tag of the PhysiCell release to use')
    group.add_argument('-pb','--pkpd_branch', default=None, help='use this branch of a PhysiCell fork rather than a release')

    parser.add_argument('--studio',action='store_true', help='also downloads a copy of studio with physipkpd integration')

def extract_from_url(url, target_dir):
    data = requests.get(url)
    local_file = append_suffix("temp",".zip")
    with open(local_file, 'wb')as f:
        f.write(data.content)

    temp_dir = append_suffix(f"{target_dir}-TEMP")
    with ZipFile(local_file, 'r') as zObject: 
        zObject.extractall(path=temp_dir)

    folder_name = os.listdir(f"./{temp_dir}")[0]
    os.rename(f"{temp_dir}/" + folder_name, target_dir)
    os.removedirs(temp_dir)
    os.remove(local_file)