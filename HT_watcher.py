"""
This is the main script, which starts calculations,
"watches" how many are running, and starts new calculations, if possible.
"""

# external imports
import os
import pickle
import numpy as np
import time
import json

# local imports
import src.utils.basic_utils as basic_utils


# number of cores per job
ncores = 16

# requested job memory (does nothing when doing local calculation)
memory = 128000  # (mb)

# select here which workflows to run for each material
script_name = ["qe_convergence_SG15", "yambo_RPA_convergence"]

# load structures, specify the dataset here
with open(os.path.join("input_data", "alexandria_filtered_4.pckl"), "rb") as f:
    data = pickle.load(f)

# if you only want to run parts of the databases, e.g. if the script crashed after some time, specify this here
import random

random.seed(0)
# random.shuffle(data)
data = data[:100]

# name of the calculation folder
calc_folder = "calc"

# get the base directory so each script knows where the source folder is
base_dir = os.getcwd()


def launch_material(material):
    # function which launches the calculations for a material

    mat_id = material.data["mat_id"]
    # start in the correct directory
    os.chdir(base_dir)

    # setup everything
    if not os.path.exists(os.path.join(calc_folder, mat_id)):
        # create directory
        os.makedirs(os.path.join(calc_folder, mat_id))

    if not os.path.exists(os.path.join(calc_folder, mat_id, "structure.pckl")):
        structure = material.structure
        name = material.data["formula"] + "_" + mat_id

        # save the structure to a file which the workflows can read
        f = open(os.path.join(calc_folder, mat_id, "structure.pckl"), "wb")
        pickle.dump([structure, name], f)
        f.close()

    if not os.path.exists(os.path.join(base_dir, "database", f"{mat_id}.json")):
        with open(os.path.join(base_dir, "database", f"{mat_id}.json"), "a") as file:
            json.dump(material.as_dict(), file)

    # create the job file and start the job
    filename = basic_utils.start_calc(
        base_dir,
        "batchjob",
        mat_id,
        script_name,
        file_name="Watcher",
        ncores=ncores,
        memory=memory,
        calc_folder=calc_folder,
        jobname="Watcher",
    )
    os.chdir(base_dir)
    return True


# loop over all materials
curr_idx = 0
while curr_idx < len(data):
    print(f"Current id: {curr_idx}")
    curr_mat = data[curr_idx]
    os.chdir(base_dir)

    # check for control instructions
    os.chdir("CONTROL")
    if os.path.isfile("STOP"):
        print("Stopping!")
        break
    if os.path.isfile("PAUSE"):
        print("Paused!")
        time.sleep(60)
        continue

    # save the current number, in case the job crashes
    np.savetxt("CURR_IDX", [curr_idx])

    # find out if any jobs can be started
    # this probably only works on our cluster and has to be adapted
    # if the script is used on another cluster
    ncores_total = np.loadtxt("NCORES")
    ncores_curr = os.popen("bjobs -sum -J Watcher").read()
    ncores_curr = ncores_curr.split("\n")
    ncores_curr = ncores_curr[1].split()
    ncores_curr = sum(int(x) for x in ncores_curr)
    max_jobs = int(ncores_total / ncores)
    curr_jobs = int(ncores_curr / ncores)
    print(curr_jobs, max_jobs)

    if curr_jobs < max_jobs:
        launch_material(curr_mat)
        curr_idx = curr_idx + 1
    else:
        print("All slots filled")
        time.sleep(60)
