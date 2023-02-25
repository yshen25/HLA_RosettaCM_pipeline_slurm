#!/usr/bin/env python3
#Shawn Aug072020
#assign most suitable template

import os
import shutil
import csv
from os import remove
from os.path import splitext
from subprocess import Popen, PIPE
from pathlib import Path
from shutil import move

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_score(seq, blastdb):
    # to build a valid blastp database, run "makeblastdb -in <template_fasta> -dbtype prot -hash_index"
    blastp = Popen(["blastp", "-db", blastdb, "-matrix", "BLOSUM90", "-outfmt", "6 sacc bitscore"], stdin=PIPE, stdout=PIPE)
    blastp_out = blastp.communicate(input=seq.encode())[0].decode().strip().split("\n")
    blastp_dict = dict([i.split("\t") for i in blastp_out])

    return blastp_dict

def template_blastp(target_seq_file, blastdb, score_file=None):

    if score_file == None:
        score_file = splitext(target_seq_file)[0] + ".score.csv"
    
    #get template ids
    templates = SeqIO.parse(blastdb, "fasta")
    template_keys = [i.id for i in templates]
    template_keys.insert(0, "Target")
    
    with open(score_file, "w") as score_handle:
        score_writer = csv.DictWriter(score_handle, template_keys)    
        score_writer.writeheader()
        for record in SeqIO.parse(target_seq_file, "fasta"):
            score_dict = get_score(str(record.seq), blastdb)
            score_dict["Target"] = record.id
            print(record.id)
            score_writer.writerow(score_dict)

    return score_file

def get_score_matrix(score_file):
    
    with open(score_file, "r", newline='') as score_handle:
        reader = csv.reader(score_handle)
        templates = next(reader)[1:]
        score_list = [np.array([int(float(j)) for j in i[1:]]) for i in reader]

    return templates, score_list

def assign_template1(target_list, score_file):
    #used for MHC I

    templates, score_list = get_score_matrix(score_file)

    templates = [i.split("_")[0] for i in templates] #remove _A chain identifier

    with open(target_list, "r", newline='') as inf, open(f"{target_list}.temp", "w", newline='') as ouf:
        reader = csv.reader(inf)
        writer = csv.writer(ouf)
        header = next(reader)
        if not "Best_template" in header:
            writer.writerow(header + ["Best_template", "Total_bitscore"])
        else:
            print("templates already assigned!")
            remove(f"{target_list}.temp")
            return 0

        for score in score_list:
            Total_bitscore = np.amax(score)
            Best_template = templates[np.argmax(score)]

            # pick 2nd best template, for validation only
            #Total_bitscore = np.sort(score)[-2]
            #Best_template = templates[np.argsort(score)[-2]]

            writer.writerow(next(reader) + [Best_template, Total_bitscore])

    remove(target_list)
    move(f"{target_list}.temp", target_list)

    return 0



def assign_template2(target_list, A_score_file, B_score_file):
    #used for MHC II

    A_templates, A_score_list = get_score_matrix(A_score_file)
    _, B_score_list = get_score_matrix(B_score_file)

    templates = [i.split("_")[0] for i in A_templates]

    with open(target_list, "r", newline='') as inf, open(f"{target_list}.temp", "w", newline='') as ouf:
        reader = csv.reader(inf)
        writer = csv.writer(ouf)
        header = next(reader)
        if not "Best_template" in header:
            writer.writerow(header + ["Best_template", "Total_bitscore"])
        else:
            print("templates already assigned!")
            remove(f"{target_list}.temp")
            return 0

        for A_score in A_score_list:
            for B_score in B_score_list:
                Total_bitscore_array = A_score + B_score
                Total_bitscore = np.amax(Total_bitscore_array)
                Best_template = templates[np.argmax(Total_bitscore_array)]

                writer.writerow(next(reader) + [Best_template, Total_bitscore])

    remove(target_list)
    move(f"{target_list}.temp", target_list)

    return 0

def blastp(query_seq, subject_file):

    blastp = Popen(["blastp", "-subject", subject_file, "-matrix", "BLOSUM90", "-outfmt", "6 qstart qend sstart qseq sseq pident"], stdin=PIPE, stdout=PIPE)
    blastp_out = blastp.communicate(input=query_seq.encode())[0].decode().strip().split("\n")[0].split("\t")
    #print(blastp_out)
    return blastp_out

def prepare_files1(target_list, seq_file, workdir=None, template_dir=None):
    #MHC I only, generate files and directory trees for building models

    if workdir == None:
        workdir = os.path.dirname(os.path.abspath(target_list))

    if template_dir == None:
        template_dir = "/home/ys0/software_local/rosetta_CM/database/Class_I/templates"
        #template_dir = "/home/shawn/work_bench/development/rosetta_CM/database/Class_I/templates"

    seq_all = list(SeqIO.parse(seq_file, "fasta"))
    with open(target_list, "r", newline='') as inf, open(f"{target_list}.temp", "w", newline='') as ouf:
        reader = csv.DictReader(inf)
        writer = csv.DictWriter(ouf, reader.fieldnames+["A_cov", "pident"])
        writer.writeheader()
        #Stat	Allele	A_index A_length    Best_template	Total_bitscore  A_cov	pident
        for row in reader:
            if not row['Stat'] == '#':
                A_gene = row['Allele']
                A_index = int(row['A_index'])-1
                template_id = row['Best_template']

                allele_dir = f"{workdir}/{A_gene}"
                
                if not os.path.exists(allele_dir):
                    os.makedirs(allele_dir)
                else:
                    print(f"{A_gene} exists!")
                    writer.writerow(row)
                    continue

                A_rec = seq_all[A_index]

                A_qstart, A_qend, A_sstart, A_qseq, A_sseq, pident = blastp(str(A_rec.seq), f"{template_dir}/{template_id}.fasta")
                # known bug in this step: gap in the begining is ignored
                # bug fixed by padding in the begining

                row["A_cov"] = f"{A_qstart}-{A_qend}"
                row["pident"] = pident

                #Full sequence
                full_rec = SeqRecord(Seq(A_qseq).ungap("-"), id=f"{A_gene}", description=A_gene)
                SeqIO.write(full_rec, f"{allele_dir}/target.fasta", "fasta")

                #Alignment files
                A_srec = SeqRecord(Seq('X' * (int(A_sstart) - 1) + A_sseq), id="template", description="template") # append X in the begining if sstart not 1
                A_qrec = SeqRecord(Seq('-' * (int(A_sstart) - 1) + A_qseq), id="target", description="target") # add gap in the begining
                SeqIO.write([A_srec, A_qrec], f"{allele_dir}/A_chain.aln", "fasta")
                
                row['Stat'] = "-"
                print(A_gene)

            writer.writerow(row)

    os.remove(target_list)
    shutil.move(f"{target_list}.temp", target_list)

    return 0

def prepare_files2(target_list, A_seq_file, B_seq_file, workdir=None, template_dir=None):
    #generate files and directory trees for building models

    if workdir == None:
        workdir = os.path.dirname(os.path.abspath(target_list))

    if template_dir == None:
        template_dir = "/home/ys0/software_local/rosetta_CM/database/Class_II/templates"

    A_seq_all = list(SeqIO.parse(A_seq_file, "fasta"))
    B_seq_all = list(SeqIO.parse(B_seq_file, "fasta"))

    with open(target_list, "r", newline='') as inf, open(f"{target_list}.temp", "w", newline='') as ouf:
        reader = csv.DictReader(inf)
        writer = csv.DictWriter(ouf, reader.fieldnames+["A_cov", "B_cov"])
        writer.writeheader()
        #Stat	Allele	A_gene	A_index	A_length	B_gene	B_index	B_length	Best_template	Total_bitscore  A_cov   B_cov
        for row in reader:
            if not row['Stat'] == '#':
                A_gene = row['A_gene']
                A_index = int(row['A_index'])-1
                B_gene = row['B_gene']
                B_index = int(row['B_index'])-1
                template_id = row['Best_template']

                allele_dir = f"{workdir}/{A_gene}/{B_gene}"
                
                if not os.path.exists(allele_dir):
                    os.makedirs(allele_dir)
                else:
                    print(f"{A_gene}_{B_gene} exists!")
                    writer.writerow(row)
                    continue

                A_rec = A_seq_all[A_index]
                B_rec = B_seq_all[B_index]

                A_qstart, A_qend, A_sstart, A_qseq, A_sseq = blastp(str(A_rec.seq), f"{template_dir}/{template_id}.fasta")
                B_qstart, B_qend, B_sstart, B_qseq, B_sseq = blastp(str(B_rec.seq), f"{template_dir}/{template_id}.fasta")

                row["A_cov"] = f"{A_qstart}-{A_qend}"
                row["B_cov"] = f"{B_qstart}-{B_qend}"

                #Full sequence of A and B chain, linked by /
                full_rec = SeqRecord(Seq(A_qseq + '/' + B_qseq).ungap("-"), id=f"{A_gene}_{B_gene}", description=f"{A_gene}_{B_gene}")
                SeqIO.write(full_rec, f"{allele_dir}/target.fasta", "fasta")
                #Full sequence of A and B chain, linked directly
                link_rec = SeqRecord(Seq(A_qseq + B_qseq).ungap("-"), id=f"{A_gene}_{B_gene}", description=f"{A_gene}_{B_gene}")
                SeqIO.write(link_rec, f"{allele_dir}/target.link.fasta", "fasta")

                #Alignment files
                A_srec = SeqRecord(Seq('X' * (int(A_sstart) - 1) + A_sseq), id="template", description="template")
                A_qrec = SeqRecord(Seq('-' * (int(A_sstart) - 1) + A_qseq), id="target", description="target")
                SeqIO.write([A_srec, A_qrec], f"{allele_dir}/A_chain.aln", "fasta")

                B_srec = SeqRecord(Seq('X' * (int(B_sstart) - 1) + B_sseq), id="template", description="template")
                B_qrec = SeqRecord(Seq('-' * (int(B_sstart) - 1) + B_qseq), id="target", description="target")
                SeqIO.write([B_srec, B_qrec], f"{allele_dir}/B_chain.aln", "fasta")
                
                row['Stat'] = "-"
                print(f"{A_gene}_{B_gene}")

            writer.writerow(row)

    os.remove(target_list)
    shutil.move(f"{target_list}.temp", target_list)

    return 0

def Class1(args):
    blastdb = str(Path(__file__).parent.parent.resolve()) + "/database/Class_I/Blast_DB/templates.fasta"
    score_file = template_blastp(args.seq_file, blastdb)
    assign_template1(args.job_file, score_file)
    prepare_files1(args.job_file, args.seq_file)
    return 0

def Class2(args):
    blastdb_A = str(Path(__file__).parent.parent.resolve()) + "/database/Class_II/Blast_DB/templates_A.fasta"
    blastdb_B = str(Path(__file__).parent.parent.resolve()) + "/database/Class_II/Blast_DB/templates_B.fasta"
    score_file_A = template_blastp(args.seq_file_A, blastdb_A)
    score_file_B = template_blastp(args.seq_file_B, blastdb_B)

    assign_template2(args.job_file, score_file_A, score_file_B)
    prepare_files2(args.job_file, args.seq_file_A, args.seq_file_B)
    return 0

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(help="choose MHC class")
    parser.set_defaults(func=lambda args: parser.print_help())

    # parser for MHC class I
    parser1 = subparser.add_parser("1", help="Class 1")
    parser1.set_defaults(func=Class1)
    parser1.add_argument("seq_file", help="sequence file generated by IMGT_input.py")
    parser1.add_argument("job_file", help="job file generated by database_access.py initial")

    # parser for MHC class II
    parser2 = subparser.add_parser("2", help="Class 2")
    parser2.set_defaults(func=Class2)
    parser2.add_argument("seq_file_A", help="sequence file for A chain")
    parser2.add_argument("seq_file_B", help="sequence file for B chain")
    parser2.add_argument("job_file", help="job file generated by database_access.py initial")

    args = parser.parse_args()
    args.func(args)
