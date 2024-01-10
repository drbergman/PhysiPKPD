import os
import argparse
import requests
from zipfile import ZipFile 
import shutil

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
    print(f"  Extracting to {local_file}...")
    with open(local_file, 'wb')as f:
        f.write(data.content)

    temp_dir = append_suffix(f"{target_dir}-TEMP")
    print(f"  Extracting to {temp_dir}...")
    with ZipFile(local_file, 'r') as zObject: 
        zObject.extractall(path=temp_dir)

    print(f"  Moving to {target_dir}...")
    folder_name = os.listdir(f"./{temp_dir}")[0]
    os.rename(f"{temp_dir}/" + folder_name, target_dir)
    os.removedirs(temp_dir)
    os.remove(local_file)

def get_studio_pkpd(USE_STUDIO, studio_dir):
    if USE_STUDIO is False:
        return
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
        extract_from_url(remote_url, studio_dir)

def get_pkpd(args, target_dir):
    if os.path.exists(f"{target_dir}/addons/PhysiPKPD"):
        print(f"Found PhysiPKPD files already in {target_dir}.")
        return
    if os.path.exists(f"{target_dir}/addons") is False:
        print(f"Making an addons folder in the project directory {target_dir}...")
        os.mkdir(f"{target_dir}/addons")

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
    else:
        remote_url = 'https://github.com/drbergman/PhysiPKPD/archive/refs/heads/' + PKPD_BRANCH + '.zip'

    print(f"PKPD remote_url = {remote_url}")
    pkpd_target_dir = append_suffix("PhysiPKPD-TEMP")
    extract_from_url(remote_url, pkpd_target_dir)
    print(f"extracted PhysiPKPD to {pkpd_target_dir}")

    os.rename(f"{pkpd_target_dir}/addons/PhysiPKPD",f"{target_dir}/addons/PhysiPKPD")
    os.rename(f"{pkpd_target_dir}/sample_projects_physipkpd",f"{target_dir}/sample_projects_physipkpd")
    os.rename(f"{pkpd_target_dir}/LICENSE",f"{target_dir}/addons/PhysiPKPD/LICENSE")
    os.rename(f"{pkpd_target_dir}/README.md",f"{target_dir}/addons/PhysiPKPD/README.md")
    shutil.rmtree(pkpd_target_dir)

    print(f"Moved PhysiPKPD files to {target_dir}. Deleted {pkpd_target_dir}")

def print_generic_pkpd_advice(i):
    print(f"\t{i}. Make sure to do the following:")
    i+=1
    print(f"\t\t* Activate Dirichlet conditions for any PK substrate (that's how they enter the microenvironnent)")
    print(f"\t\t* Set nonzero uptake rates for every (cell type, substrate) pairing that has a PD model (that's how damage is accumulated)")
    print(f"\t\t* Add rules for how `custom:S_damage` (S = substrate name) affects target cell types; enable rules")
    print(f"\t{i}. Consider using `damage_coloring` to assess your PD model parameters as you build the model.")
    i+=1
    print(f"\t     In custom_modules/custom.cpp, replace my_coloring_function with\n")
    print("std::vector<std::string> my_coloring_function(Cell *pC)\n{ return damage_coloring(pC); }\n\n")


def update_physicell_files(DIR, project_loaded=False):
    main_file = f'{DIR}/main.cpp'
    if project_loaded and os.path.exists(main_file) is False:
        print(f"{main_file} does not exist so the addition of PhysiPKPD is incomplete. Make your project and re-run this to complete.")
        exit(-1)

    with open(f'{DIR}/addons/PhysiPKPD/Makefile-PhysiPKPD-Samples.txt', "r") as f:
        add_samples_to_makefile(f, f"{DIR}/sample_projects/Makefile-default")

    if project_loaded:
        makefile = f'{DIR}/Makefile'
        makefile_temp = append_suffix(f'{DIR}/Makefile-TEMP')
        add_to_makefile = True
        with open(makefile, 'r') as f:
            lines = f.readlines()
            for line_number, line in enumerate(lines):
                if "PhysiPKPD_OBJECTS := PhysiPKPD_PK.o PhysiPKPD_PD.o" in line:
                    print(f"WARNING: Found PhysiPKPD_OBJECTS already defined in {makefile}. Skipping the rest of Makefile additions.")
                    print(f"\tIf compilation errors occur, make sure PhysiPKPD_PK.o and PhysiPKPD_PD.o are defined")
                    print(f"\tAnd make sure that $(PhysiPKPD_OBJECTS) is added to PhysiCell_OBJECTS.")
                    add_to_makefile = False
                    break
            if add_to_makefile:
                for line_number, line in enumerate(lines):
                    if "PhysiCell_custom_module_OBJECTS :=" in line:
                        lines.insert(line_number,"PhysiPKPD_OBJECTS := PhysiPKPD_PK.o PhysiPKPD_PD.o\n\n")
                        break
        
        if add_to_makefile:
            source_file = f'{DIR}/addons/PhysiPKPD/Makefile-PhysiPKPD-Objects.txt'
            with open(source_file, 'r') as sf:
                new_lines = sf.readlines()

            for line_number, line in enumerate(lines):
                if "custom_modules/custom.cpp" in line:
                    break

            with open(makefile_temp,'w') as f:
                for i in range(line_number):
                    print(lines[i].rstrip('\n'), file=f)
                for new_line in new_lines:
                    print(new_line.rstrip('\n'), file=f)
                print("\n", file=f)
                for i in range(i+1,len(lines)):
                    print(lines[i].rstrip('\n'), file=f)

            with open(makefile_temp,'r+') as f:
                lines = f.readlines()

            for i, line in enumerate(lines):
                if "PhysiCell_OBJECTS :=" in line:
                    line = line.rstrip('\n') + ' $(PhysiPKPD_OBJECTS)\n'
                    break
            lines[i] = line

            with open(makefile_temp,'w') as f:
                f.writelines(lines)

            print(f"{makefile_temp} is ready for PhysiPKPD. Editing other files first before overwriting {makefile}")

        def get_line_number(s, lines, print_missing_message=True):
            for line_number, line in enumerate(lines):
                if s in line:
                    return line_number
            if print_missing_message:
                print(f"{f.name} does not have a line containing {s}")
            return None
            
        add_to_main = True
        main_temp = append_suffix(f'{DIR}/main','.cpp')
        with open(main_file, 'r+') as f:
            lines = f.readlines()
            if get_line_number("setup_pharmacodynamics", lines, print_missing_message=False) is not None:
                print(f"WARNING: Found setup_pharmacodynamics already defined in {main_file}. Skipping the rest of main.cpp additions.")
                print(f"\tIf errors occur, make sure PK_model(PhysiCell_globals.current_time); and PD_model(PhysiCell_globals.current_time); are on either side of microenvironment.simulate_diffusion_decay.")
                add_to_main = False
            else:
                line_number = get_line_number("setup_tissue", lines)
                if line_number is not None:
                    lines.insert(line_number+1, "\n\tsetup_pharmacodynamics();\n")  # new_string should end in a newline
                    with open(main_temp,'w') as f_temp:
                        f_temp.writelines(lines)  # No need to truncate as we are increasing filesize
                else:
                    print(f"{DIR}/main.cpp does not include `setup_tissue();`???")
                    exit()
                
        if add_to_main:
            with open(main_temp, 'r+') as f:
                lines = f.readlines()
                line_number = get_line_number("microenvironment.simulate_diffusion_decay", lines)
                if line_number is not None:
                    lines.insert(line_number+1, "\n\t\t\tPD_model( PhysiCell_globals.current_time );\n")  # new_string should end in a newline
                    lines.insert(line_number, "\t\t\tPK_model( PhysiCell_globals.current_time );\n\n")  # new_string should end in a newline
                    f.seek(0)  # readlines consumes the iterator, so we need to start over
                    f.writelines(lines)  # No need to truncate as we are increasing filesize
                else:
                    print(f"{DIR}/main.cpp does not include `microenvironment.simulate_diffusion_decay(diffusion_dt);`???")
                    exit()

            print(f"{main_temp} is ready for PhysiPKPD. Editing other files first before overwriting {main_file}")

        with open(f"{DIR}/custom_modules/custom.h", "r+") as f:
            lines = f.readlines()
            if get_line_number("addons/PhysiPKPD/src/PhysiPKPD.h", lines, print_missing_message=False):
                print("Found addons/PhysiPKPD/src/PhysiPKPD.h in {DIR}/custom_modules/custom.h. Not going to add it again here.")
            else:
                line_number = get_line_number("#include", lines)
                if line_number is not None:
                    lines.insert(line_number,'#include "../addons/PhysiPKPD/src/PhysiPKPD.h"\n')
                    f.seek(0)
                    f.writelines(lines)
                else:
                    print(f"{DIR}/custom_modules/custom.h does not have any `#include` statements???")
                    exit()

        if add_to_makefile:
            with open(makefile, 'w') as f:
                with open(makefile_temp, 'r') as f_temp:
                    lines = f_temp.readlines()
                    f.writelines(lines)

        if add_to_main:
            with open(main_file, 'w') as f:
                with open(main_temp, 'r') as f_temp:
                    lines = f_temp.readlines()
                    f.writelines(lines)

        print(f"All files (Makefile, main.cpp, and custom_modules/custom.h) updated and overwritten now. Removing temporary files...")
        if os.path.exists(makefile_temp):
            os.remove(makefile_temp)
        if os.path.exists(main_temp):
            os.remove(main_temp)
    else:
        with open(f'{DIR}/addons/PhysiPKPD/Makefile-PhysiPKPD-Samples.txt', "r") as f:
            add_samples_to_makefile(f, f"{DIR}/Makefile")

def add_samples_to_makefile(samples_file, path_to_makefile):
    # Update Makefile
    print("----------------------")
    print(f"Now updating {path_to_makefile} to be ready to make PhysiPKPD samples and projects.")

    with open(path_to_makefile, 'a') as f:
        f.write("\n")
        shutil.copyfileobj(samples_file, f)
