#!/usr/bin/env python3
"""
    read the project file, and extract top models
"""
# Shawn 112219

import argparse
import os
import shutil

from database_access import StatusEditor, read_csv, allele_lookup
from pick_model import extract, pbs_relax


class picker():
    def __init__(self) -> None:
        self.extracts = []
        pass

    def extract_from_record(self, row, output="output", scorefile="scorefile"):
        # extract one allele from one row of record list
        TargetAllele = row['Allele']

        full_path = f"{self.work_dir}/{TargetAllele}/models"
        
        extracted_tags = extract(f"{full_path}/{output}", f"{full_path}/{scorefile}", job_prefix=TargetAllele)
        extracted_pdbs = [TargetAllele + s +".pdb" for s in extracted_tags]
        self.extracts += extracted_pdbs
        
        return

    def move_pdbs(self):
        for file in self.extracts:
            shutil.move(file, self.out_dir)
            print(file)
        print("====================================")
        print(f"extract destination: {self.out_dir}")
        return

    def pick(self, JobFile, TargetList=None, work_dir=None, out_dir=None):
        # HLA class 1 only
        if work_dir is None:
            self.work_dir = os.path.dirname(JobFile)
        else:
            self.work_dir = work_dir

        if out_dir is None:
            self.out_dir = "finished_pdbs"
        else:
            self.out_dir = out_dir

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        record_list = StatusEditor(JobFile, work_dir).check_stat(TargetList, return_list=1)

        for row in record_list:
            if row['Stat'] == "+":
                self.extract_from_record(row)

        self.move_pdbs()

        return

def direct_pick(args):
    picker().pick(args.project_file, args.list, args.work_dir, args.out_dir)
    return

def relax(args):
    # undone
    pbs_relax(args.work_dir)
    return

def relax_pick(args):
    # undone
    picker().pick(args.project_file, args.list, args.work_dir, args.out_dir)
    return

if __name__ == "__main__":

    # parser = argparse.ArgumentParser()
    # parser.set_defaults(func=lambda args: parser.print_help())
    # subparser = parser.add_subparsers(help="choose function")

    # # parser for direct_relax command
    # parser_a = subparser.add_parser("direct_pick", help="pick top model(s) from output, without relaxation)")
    # parser_a.set_defaults(func=direct_pick)
    # parser_a.add_argument("project_file", help="project csv file name")
    # parser_a.add_argument('-l', '--list', help="list file that contains target alleles, one allele each line. if not provided, extract all")
    # parser_a.add_argument('--work_dir', help="the root directory of project, default is where the project file exist")
    # parser_a.add_argument('--out_dir', help="directory that the extracted pdb files will be copied to")
    
    # # parser for relax function
    # parser_b = subparser.add_parser("relax", help="relax top models in output using FastRelax")
    # parser_b.set_defaults(func=relax)
    # parser_b.add_argument("project_file", help="project csv file name")
    # parser_a.add_argument('-l', '--list', help="list file that contains target alleles, one allele each line. if not provided, extract all")
    # parser_b.add_argument('-n', '--num', help="number of relax jobs to run")

    # # parser for relax_pick function
    # parser_c = subparser.add_parser("relax_pick", help="pick top model(s) from relax.out, after relaxation")
    # parser_c.set_defaults(func=relax_pick)
    # parser_c.add_argument("project_file", help="project csv file name")
    # parser_c.add_argument('-l', '--list', help="list file that contains target alleles, one allele each line. if not provided, extract all")
    # parser_c.add_argument('--work_dir', help="the root directory of project, default is where the project file exist")
    # parser_c.add_argument('--out_dir', help="directory that the extracted pdb files will be copied to")

    # args = parser.parse_args()
    # args.func(args)

    parser = argparse.ArgumentParser()
    parser.add_argument('project_file', help= "project file name")
    parser.add_argument('-l', '--list', help="list file that contains target alleles, one allele each line. if not provided, extract all")
    parser.add_argument('--work_dir', help="the root directory of project, default is where the project file exist")
    parser.add_argument('--out_dir', help="directory that the extracted pdb files will be copied to")
    args = parser.parse_args()
    picker().pick(args.project_file, args.list, args.work_dir, args.out_dir)