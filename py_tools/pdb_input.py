#!/usr/bin/env python3
#Shawn Jul152020
#processes structure from IMGT/3Dstructure-DB
import os
import glob

from Bio import SeqIO

import extract_chain

def IMGT_filename(PDB_id):
    filename = f"IMGT-{PDB_id.upper()}.pdb"
    return filename

def template_info_read(template_list):
    #read cleaned template pdb from csv file
    import csv
    
    template_pdbs = []
    with open(template_list, newline='') as csv_handle:
        template_reader = csv.reader(csv_handle)
        next(template_reader)
        for row in template_reader:
            if not row[0].startswith("#"):
                template_pdbs.append(row[0])
    return template_pdbs

def combine_template_seq(template_D, template_list, output_prefix="templates"):

    template_ids = template_info_read(template_list)

    A_collection = []
    B_collection = []

    for pid in template_ids:
        record = SeqIO.parse(f"{template_D}/{pid}.fasta", "fasta")
        A_collection.append(next(record))
        B_collection.append(next(record))

    with open(f"{output_prefix}_A.fasta", "w") as ouf:
        SeqIO.write(A_collection, ouf, "fasta")
    print(f"{output_prefix}_A.fasta")

    with open(f"{output_prefix}_B.fasta", "w") as ouf:
        SeqIO.write(B_collection, ouf, "fasta")
    print(f"{output_prefix}_B.fasta")
    
    return 0
"""
def main(raw_D, template_D, template_list):

    template_pdbs = filter(None, template_info_read(template_list))

    for pid in template_pdbs:
        pdb_infile = f"{raw_D}/{IMGT_filename(pid)}"
        if not os.path.exists(pdb_infile):
            raise ValueError(f"PDB file {IMGT_filename(pid)} not exist!!")
        else:
            print(pid)
            out_filename = f"{template_D}/{pid}.pdb"
            extract_chain.extract(pdb_infile, "AB", out_filename)
            extract_chain.extract_sequence(out_filename)
            
    combine_template_seq(template_D, template_list)

    return 0
"""

def main(raw_D, out_D, MHC_class):
    """
    clean pdb files downloaded from IMGT-3D database
    raw_D: directory that contains downloaded pdb files
    out_D: directory that is used to store cleaned files
    """

    if not os.path.exists(out_D):
        os.makedirs(out_D)

    for in_file in glob.glob(f"{raw_D}/*.pdb"):
        out_filename = f"{out_D}/{os.path.basename(in_file)}"
        if MHC_class == "1":
            extract_chain.extract(in_file, "A", out_filename)
            extract_chain.extract_sequence(out_filename)

        elif MHC_class == "2":
            extract_chain.extract(in_file, "AB", out_filename)
            extract_chain.extract_sequence(out_filename)
        print(f"processed: {in_file}")

    return

    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="clean pdb files from IMGT-3D structure to use as template")
    parser.add_argument("input_dir", help="directory containing downloaded pdb files")
    parser.add_argument("output_dir", help="directory for cleaned output files")
    parser.add_argument("MHC_class", choices=["1","2"], help="MHC class, 1 or 2")
    #parser.add_argument("template_list", help="csv file with info from IMGT-3D structure")
    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.MHC_class)