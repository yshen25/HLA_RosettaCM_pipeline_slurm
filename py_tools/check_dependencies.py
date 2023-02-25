#!usr/bin/env python3
#check dependencies
#shawn 081419
from shutil import which
from os import path

def check_dependencies(template_db, template_path, psipred_path, psiblast_db, fragwghts, wts1, wts2, wts3, vall, psiblast, blastp, mafft):
    status = 0
    if which("blastp") is None and not path.exists(blastp):
        print("BLAST+ package not found!")
        status = 1
    else:
        print("...BLAST+ PASS")
    if which("mafft") is None and not path.exists(mafft):
        print("mafft not found!")
        status = 1
    else:
        print("...mafft PASS")
    if which("psiblast") is None and not path.exists(psiblast):
        print("psiblast not found!")
        status = 1
    else:
        print("...psiblast PASS")
    if not path.exists(template_db):
        print("template_db not exist!")
        status = 1
    if not path.exists(template_path):
        print("template_path not exist!")
        status = 1
    if not path.exists(psipred_path):
        print("psipred_path not exist!")
        status = 1
    if not path.exists(fragwghts):
        print("fragment weight file not exist!")
        status = 1
    if not path.exists(wts1) or not path.exists(wts2) or not path.exists(wts3):
        print("hybridize mover weight file not exist!")
        status = 1
    if not path.exists(vall):
        print("vall database not exist!")
        status = 1
    #final status
    if status == 1:
        raise Exception('Denpendencies not satisfied! Please check PATH')
    else:
        print("ALL PASS")
    return 0
