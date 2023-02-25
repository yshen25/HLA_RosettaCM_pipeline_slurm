#!/usr/bin/env python3
from subprocess import Popen, PIPE, DEVNULL
import os
import tempfile
import sys
from pathlib import Path
from os.path import splitext, dirname
import chkparse

def MSApred(target_fasta, msa_file, msa_target_idx, coverage, psipred_path, work_dir=None):
    #psipred
    #modified from original psipred bash script
    #target_fasta is not actually used, but fed psiblast for convinience
    if work_dir == None:
        work_dir = os.path.dirname(os.path.abspath(target_fasta))

    pssm = tempfile.NamedTemporaryFile()
    #------------experimental feature---------------
    #use MSA instead of psiblast to accelerate searching
    PSI = Popen(["psiblast", "-subject", target_fasta, "-in_msa", msa_file, "-msa_master_idx", msa_target_idx, "-out_pssm", pssm.name], stdout=DEVNULL) 
    PSI.communicate()

    mtx = tempfile.NamedTemporaryFile()

    checkpoint = f"{work_dir}/frag.checkpoint" # final output 1
    chkparse.single(pssm.name, mtx.name, checkpoint, coverage) # coverage is 1-based, and stop index is included

    ss = tempfile.NamedTemporaryFile()
    psipred = Popen([f"{psipred_path}/bin/psipred", mtx.name, f"{psipred_path}/data/weights.dat", f"{psipred_path}/data/weights.dat2", f"{psipred_path}/data/weights.dat3"], stdout=PIPE)
    with open(ss.name, "w") as sshandle:
        sshandle.write(psipred.communicate()[0].decode())

    ss2 = f"{work_dir}/frag.ss2" # final output 2
    psipass = Popen([f"{psipred_path}/bin/psipass2", f"{psipred_path}/data/weights_p2.dat" ,"1" ,"1.0" ,"1.0" ,ss2, ss.name], stdout=DEVNULL)
    psipass.communicate()

    return checkpoint, ss2

def MSApred_multichain(target_fasta, msa_file_list:list, msa_target_idx_list:list, coverage_list:list, psipred_path, work_dir = None):
    #psipred
    #modified from original psipred bash script, and altered to process multi-chain protein
    if work_dir == None:
        work_dir = os.path.dirname(os.path.abspath(target_fasta))
    
    #one combined midterm result
    mtx = tempfile.NamedTemporaryFile()
    #two final output files
    checkpoint = f"{work_dir}/frag.checkpoint"
    ss2 = f"{work_dir}/frag.ss2"
    pssm = []

    for msa_file, msa_target_idx in zip(msa_file_list, msa_target_idx_list):
        pssm.append(tempfile.NamedTemporaryFile().name)
        PSI = Popen(["psiblast", "-subject", target_fasta, "-in_msa", msa_file, "-msa_master_idx", msa_target_idx, "-out_pssm", pssm[-1]], stdout=DEVNULL, stderr=DEVNULL) 
        PSI.communicate()

    chkparse.double(pssm[0], pssm[1], mtx.name, checkpoint, coverage_list[0], coverage_list[1])

    #the pipeline runs as usual
    ss = tempfile.NamedTemporaryFile()
    psipred = Popen([f"{psipred_path}/bin/psipred", mtx.name, f"{psipred_path}/data/weights.dat", f"{psipred_path}/data/weights.dat2", f"{psipred_path}/data/weights.dat3"], stdout=PIPE)
    with open(ss.name, "w") as sshandle:
        sshandle.write(psipred.communicate()[0].decode())

    ss2 = f"{work_dir}/frag.ss2"
    psipass = Popen([f"{psipred_path}/bin/psipass2", f"{psipred_path}/data/weights_p2.dat" ,"1" ,"1.0" ,"1.0" ,ss2, ss.name], stdout=DEVNULL)
    psipass.communicate()
    
    return checkpoint, ss2

def frag_flag(target_fasta, checkpoint, ss2, fragwghts, vall):
    #write flags for frag_picker
    work_dir = os.path.dirname(os.path.abspath(target_fasta))
    flags = f"{work_dir}/frag.flag"
    #check if all needed files exists
    if not os.path.exists(target_fasta):
        sys.exit("target fasta file not exists")
    if not os.path.exists(checkpoint):
        sys.exit("checkpoint file not exists")
    if not os.path.exists(ss2):
        sys.exit("ss2 file not exists")
    if not os.path.exists(fragwghts):
        sys.exit("weights file not exists")
    if not os.path.exists(work_dir):
        sys.exit("working path not created")
    
    #write the flag
    with open(flags, "w") as fh:
        fh.write(f"-in::file::vall\t{vall}\n")
        fh.write(f"-in::file::checkpoint\t{checkpoint}\n")
        fh.write(f"-in::file::fasta\t{target_fasta}\n")
        fh.write(f"-frags::ss_pred\t{ss2} psipred\n")
        fh.write(f"-frags::scoring::config\t{fragwghts}\n")
        fh.write(f"-frags::bounded_protocol\n")
        fh.write("-frags::frag_sizes\t3 9\n-frags::n_candidates\t200\n-frags::n_frags\t200\n")
        fh.write(f"-out::file::frag_prefix\t{work_dir}/target\n")
        fh.write(f"-frags::describe_fragments\t{work_dir}/target.fsc\n")
        fh.write("-overwrite")
    return flags
#run frag_picker

def pbs_run(target_fasta, msa_file, msa_target_idx, coverage, psiblast_db = None, multi_chain = 0):
    #generate other names and paths
    #msa index is 1 based, used by psiblast
    #coverage is 1 based, generated by blastp
    if psiblast_db == None:
        psiblast_db = "/home/ys0/software_local/lib/uniref50/uniref50.fasta"
    psipred_path = "/home/ys0/software_local/psipred_4.02"
    fragwghts = "/home/ys0/software_local/rosetta_CM/parameters/FragPick.wghts"
    vall = "/home/ys0/software_local/rosetta_2019.22.60749/tools/fragment_tools/vall.jul19.2011"
    pbs_dir = str(Path(__file__).parent.parent.resolve())+"/pbs_scripts"
    work_dir = dirname(target_fasta)

    if multi_chain:
        checkpoint, ss2 = MSApred_multichain(target_fasta, msa_file, msa_target_idx, coverage, psipred_path)
    else:
        checkpoint, ss2 = MSApred(target_fasta, msa_file, msa_target_idx, coverage, psipred_path)
    
    flag_file = frag_flag(target_fasta, checkpoint, ss2, fragwghts, vall)

    pbs_submit = Popen(["sbatch", f"{pbs_dir}/2.frag_picker.pbs", flag_file, work_dir], stdout=PIPE, stderr=PIPE)
    ticket, error = [x.decode() for x in pbs_submit.communicate()]
    if error:
        print(error)
        status = 1
    else:
        status = 0
    return ticket, status

def direct_run(target_fasta, psiblast_db = "/home/ys0/software_local/lib/uniref50/uniref50.fasta"):
    #used to run fragment picker interactively
    psipred_path = "/home/ys0/software_local/psipred_4.02"
    fragwghts = "/home/ys0/software_local/rosetta_CM/parameters/FragPick.wghts"
    vall = "/home/ys0/software_local/rosetta_2019.22.60749/tools/fragment_tools/vall.jul19.2011"
    pbs_dir = str(Path(__file__).parent.parent.resolve())+"/pbs_scripts"

    chk = tempfile.NamedTemporaryFile()
    print("running psiblast...")
    PSI = Popen(["psiblast", "-num_threads", "4", "-db", psiblast_db, "-query", target_fasta, "-inclusion_ethresh", "0.0001", "-out_pssm", chk.name, "-num_iterations", "3", "-num_alignments", "0", "-num_descriptions", "500"], stdout=DEVNULL, stderr=DEVNULL)
    PSI.communicate()
    print("...psiblast done")
    mtx = tempfile.NamedTemporaryFile()
    print("predicting secondary structure...")
    checkpoint = f"{splitext(target_fasta)[0]}.checkpoint"
    chkparse = Popen([f"{psipred_path}/bin/chkparse_alt", chk.name, checkpoint], stdout=PIPE)
    with open(mtx.name, "w") as mtxhandle:
        mtxhandle.write(chkparse.communicate()[0].decode())
    print("   Pass1...")
    
    ss = tempfile.NamedTemporaryFile()
    psipred = Popen([f"{psipred_path}/bin/psipred", mtx.name, f"{psipred_path}/data/weights.dat", f"{psipred_path}/data/weights.dat2", f"{psipred_path}/data/weights.dat3"], stdout=PIPE)
    with open(ss.name, "w") as sshandle:
        sshandle.write(psipred.communicate()[0].decode())
    print("   Pass2...")

    ss2 = f"{splitext(target_fasta)[0]}.ss2"
    psipass = Popen([f"{psipred_path}/bin/psipass2", f"{psipred_path}/data/weights_p2.dat" ,"1" ,"1.0" ,"1.0" ,ss2, ss.name], stdout=DEVNULL)
    psipass.communicate()
    print(f"...psipred done")
    print("picking fragments...")
    flag_file = frag_flag(target_fasta, checkpoint, ss2, fragwghts, vall)
    pbs = Popen(["sbatch","--account","bsd","--partition","testing","--mem=16G","-t","01:30:00","--nodes","1","--ntasks","16",f"{pbs_dir}/frag_picker.sh", f"{flag_file}"])
    #pbs = Popen(["/home/ys0/software_local/rosetta_CM/pbs_scripts/frag_picker.sh", f"{flag_file}"])
    pbs.communicate()
    print("...rosetta done")

    return 0

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="fragment picker")
    parser.add_argument("target_fasta")
    args = parser.parse_args()
    direct_run(args.target_fasta)