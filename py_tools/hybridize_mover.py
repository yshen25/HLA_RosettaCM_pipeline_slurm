#!/usr/bin/env python3
from subprocess import Popen, PIPE
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
import os

def hybridize_xml(work_dir, thread_pdb, frag_3mer, frag_9mer, wts_trilogy):
    xml = f"{work_dir}/hybridize.xml"
    
    with open(xml, "w") as xml_handle:
        xml_content = f'''<ROSETTASCRIPTS>
    <TASKOPERATIONS>
    </TASKOPERATIONS>
    <SCOREFXNS>
        <ScoreFunction name="stage1" weights="{wts_trilogy[0]}" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="1"/>
        </ScoreFunction>
        <ScoreFunction name="stage2" weights="{wts_trilogy[1]}" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="0.5"/>
        </ScoreFunction>
        <ScoreFunction name="fullatom" weights="{wts_trilogy[2]}" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="0.5"/>
        </ScoreFunction>
    </SCOREFXNS>
    <FILTERS>
    </FILTERS>
    <MOVERS>
        <Hybridize name="hybridize" stage1_scorefxn="stage1" stage2_scorefxn="stage2" fa_scorefxn="fullatom" batch="1" stage1_increase_cycles="1.0" stage2_increase_cycles="1.0" linmin_only="1">
            <Fragments three_mers="{frag_3mer}" nine_mers="{frag_9mer}"/>
            <Template pdb="{thread_pdb}" cst_file="AUTO" weight="1.000" />
        </Hybridize>
    </MOVERS>
    <APPLY_TO_POSE>
    </APPLY_TO_POSE>
    <PROTOCOLS>
        <Add mover="hybridize"/>
    </PROTOCOLS>

</ROSETTASCRIPTS>'''
        xml_handle.write(xml_content)
    return xml


def hybridize_flag(target_fasta, xml, work_dir):
    flag_file = f"{work_dir}/hybridize.flag"
    with open(flag_file, "w") as flag_handle:
        flag_content = f'''-in:file:fasta {target_fasta}
-parser:protocol {xml}
-nstruct 10
-out:path:all {work_dir}/models
-out:file:silent output
-out:file:scorefile scorefile
-relax:minimize_bond_angles
-relax:minimize_bond_lengths
-relax:jump_move true
-default_max_cycles 200
-relax:min_type lbfgs_armijo_nonmonotone
-score:weights ref2015_cart
-use_bicubic_interpolation
-hybridize:stage1_probability 1.0
-overwrite'''
        flag_handle.write(flag_content)
    
    return flag_file

def pbs_run(target_fasta, thread_pdb, frag_3mer, frag_9mer, frag_ticket, work_dir=None, wts_trilogy=None):
    target_fasta = os.path.realpath(target_fasta)
    thread_pdb = os.path.realpath(thread_pdb)
    frag_3mer = os.path.realpath(frag_3mer)
    frag_9mer = os.path.realpath(frag_9mer)

    if work_dir == None:
        work_dir = os.path.dirname(target_fasta)
    else:
        work_dir = os.path.realpath(work_dir)

    if wts_trilogy == None:
        para_dir = str(Path(__file__).parent.parent.resolve())+"/parameters"
        wts_trilogy = [f"{para_dir}/stage1.wts", f"{para_dir}/stage2.wts", f"{para_dir}/stage3.wts"]

    if not os.path.exists(f"{work_dir}/models"):
        os.makedirs(f"{work_dir}/models")
    
    xml = hybridize_xml(work_dir, thread_pdb, frag_3mer, frag_9mer, wts_trilogy)
    flag_file = hybridize_flag(target_fasta, xml, work_dir)

    pbs_dir = str(Path(__file__).parent.parent.resolve())+"/pbs_scripts"
    pbs = Popen(["sbatch", f"--dependency=afterok:{frag_ticket}", f"{pbs_dir}/3.hybridize_mover.pbs", flag_file, work_dir], stdout=PIPE, stderr=PIPE)
    ticket, error = [x.decode() for x in pbs.communicate()]
    if error:
        print(error)
        status = 1
    else:
        status = 0
    return ticket, status

def direct_run(target_fasta, thread_pdb, frag_3mer, frag_9mer, work_dir=None, wts_trilogy=None):
    if work_dir is None:
        work_dir = os.path.dirname(os.path.realpath(target_fasta))
    if wts_trilogy is None:
        para_dir = str(Path(__file__).parent.parent.resolve())+"/parameters"
        wts_trilogy = [f"{para_dir}/stage1.wts", f"{para_dir}/stage2.wts", f"{para_dir}/stage3.wts"]
    print(work_dir)
    if not os.path.exists(f"{work_dir}/models"):
        os.makedirs(f"{work_dir}/models")
    
    xml = hybridize_xml(work_dir, thread_pdb, frag_3mer, frag_9mer, wts_trilogy)
    flag_file = hybridize_flag(target_fasta, xml, work_dir)

    pbs_dir = str(Path(__file__).parent.parent.resolve())+"/pbs_scripts"
    pbs = Popen(["sbatch", f"--export=F={flag_file}", "-J", "hybridize_mover", f"{pbs_dir}/3.hybridize_mover.pbs"])
    pbs.communicate()

    return 0

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="hybridize mover")
    parser.add_argument("target_fasta")
    parser.add_argument("thread_pdb")
    parser.add_argument("frag_3mer")
    parser.add_argument("frag_9mer")
    args = parser.parse_args()
    direct_run(args.target_fasta, args.thread_pdb, args.frag_3mer, args.frag_9mer)
