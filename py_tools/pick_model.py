#!/usr/bin/env python3
"""
    read the silent file, and 
    relaxation is optional
    extract top models
"""
# Shawn 112519

import re, sys, os
from subprocess import Popen, PIPE
from pathlib import Path
import argparse

def pick_top(score_file, num_top):
    top_models = []
    SORT = Popen(["sort", "-g", "-k", "2", score_file], stdout = PIPE)
    sf_handle = SORT.communicate()[0].decode()
    for line in sf_handle.splitlines()[2:num_top+2]:
        tag = re.split(" +", line.strip())[-1]
        top_models.append(tag)
    return top_models

def extract(silent_file, score_file, num_top=1, job_prefix=None):
    if job_prefix is None:
        job_prefix = silent_file.split(".")[0]
    
    top_models = pick_top(score_file, num_top)
    
    #print(top_models)
    pbs_dir = str(Path(__file__).parent.parent.resolve())+"/pbs_scripts"
    extract = Popen([f"{pbs_dir}/5.extract_pdb.sh", f"{job_prefix}", f"{silent_file}", ' '.join(top_models)])
    extract.communicate()
    return top_models

# ===== not complete =====
def pbs_relax(work_dir, num_top=3, num_per_top=3):

    top_models = pick_top(f"{work_dir}/models/scorefile", num_top)
    pbs_dir = str(Path(__file__).parent.parent.resolve())+"/pbs_scripts"

    flags = f"-parser:protocol {pbs_dir}/FastRelax.xml -nstruct {num_per_top} -in:file:silent {work_dir}/models/output -in:file:tags {' '.join(top_models)} -out:file:silent relax.out -out:file:scorefile relax.sc"
    pbs = Popen(["sbatch", f"{pbs_dir}/4.relax.pbs", flags, work_dir])
    pbs.communicate()

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("silent_file")
    parser.add_argument("score_file")
    parser.add_argument("num_top", type=int, nargs='?', default=1)
    parser.add_argument("-n", "--name", type=str)
    args = parser.parse_args()
    extract(args.silent_file, args.score_file, args.num_top, args.name)