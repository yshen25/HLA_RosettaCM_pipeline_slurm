#!/usr/bin/env python3
#Shawn Jul132020
#This module processes IMGT alignment file that contains sequence and allele info
import re
import argparse
import tempfile
from Bio import AlignIO, SeqIO, Seq
from numpy import dtype


def leading_pep(IMGT_file):
    #get the length of leading peptide
    with open(IMGT_file, "r") as inf:
        for _ in range(8): #retrieve the 9th line
            inf.readline()
        mp_start = inf.readline().find("|", 21)
        leading_peptide = inf.readline()[20:mp_start]
    lpep_length = len(leading_peptide.replace(" ", ""))
    return lpep_length

def rename_allele(allele_full):
    #rename alleles
    #only keeps two fields, ignore synonymous variations
    #and change delimiters
    allele_short = (':'.join(allele_full.split(':')[0:2]))
    return allele_short

def convert(IMGT_file):
    #convert IMGT file to clustal w
    
    def extract_refseq(inf):
        #read the leading sequence as reference to decypher following sequence
        data = inf.read()
        section_pattern = re.compile(r"(?:Prot.*\n.*\n)\s\w+\*\S{2,15}\s*(.*)", re.MULTILINE)
        template = []
        for item in re.finditer(section_pattern, data):
            template.append(item.group(1).replace(' ', '').replace('.', '-'))
        inf.seek(0)
        return template
    
    clustal_temp = tempfile.NamedTemporaryFile().name
    loop_index = -1

    with open(IMGT_file, "r") as inf, open(clustal_temp, "w") as outf:
        
        template = extract_refseq(inf)
        
        outf.write('CLUSTAL W\n')
        for line in inf:
            if line.startswith('#') or line.startswith('  '): #jump headers or empty lines
                continue
            elif line.startswith(' Prot'): #start a round
                outf.write('\n')
                loop_index += 1
            elif not line.startswith(' '): #EOF
                break
            else:
                line1 = line.strip().split(' ', 1)
                #print(line1)
                if line1[0].endswith(('N', 'L', 'S', 'C', 'A', 'Q')):
                    continue #skip alleles that doesn't express
                else:
                    try:
                        old_rec = line1[1].replace(' ', '')
                    except:
                        old_rec = ''
                    final_line = ''
                    sequence = ''
                    final_line += f"{line1[0]: <18}"
                    for i in range(len(old_rec)):
                        if old_rec[i] == '-':
                            sequence += template[loop_index][i]
                        elif (old_rec[i] == '.') or (old_rec[i] == '*'):
                            sequence += '-'
                        elif old_rec[i] == 'X':
                            break
                        else:
                            sequence += old_rec[i]
                    final_line += f"{sequence:{'-'}<{len(template[loop_index])}}"
                    final_line += '\n'
                    outf.write(final_line)
    return clustal_temp

def record_info(seq_file):
    #conclude allele info into csv table
    import csv

    out_csv = seq_file.rsplit('.', 1)[0] + '.csv'
    sequences = SeqIO.parse(seq_file, "fasta")

    with open(out_csv, "w", newline='') as outf:
        info_writer = csv.writer(outf, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        info_writer.writerow(["Allele_name", "Peptide_length", "Sequence_hash"])
        for record in sequences:
            Allele_name = record.id
            Peptide_length = len(record)
            Sequence_hash = hash(record.seq)
            info_writer.writerow([Allele_name, Peptide_length, Sequence_hash])
    print(f"Allele info concluded in: {out_csv}")
    return 0

def main(IMGT_file, out_pre = None, cut_tail_position:int=0):
    """
    cut tail position is the 1-based start index of tail, which is not included in output
    """
    if out_pre is None:
        out_pre = IMGT_file.split('.')[0] + '.nr'

    clustal_inf = convert(IMGT_file)

    lpep_length = leading_pep(IMGT_file)
    
    seq_outf = out_pre + ".faa"
    aln_outf = out_pre + ".aln"

    with open(clustal_inf, "r") as inf, open(seq_outf, "w") as soutf, open(aln_outf, "w") as aoutf:
        alignments = AlignIO.read(inf, 'clustal')
        HLA_fam_appeared = []
        for record in alignments:
            HLA_fam = rename_allele(record.id)
            #HLA_fam = record.id
            if not HLA_fam in HLA_fam_appeared:
                record.id = HLA_fam
                record.seq = record.seq[lpep_length:] #cut the leading peptide
                if cut_tail_position:
                    record.seq = record.seq[:cut_tail_position-1]
                SeqIO.write(record, aoutf, 'fasta')
                record.seq = record.seq.ungap('-')
                SeqIO.write(record, soutf, 'fasta')
                HLA_fam_appeared.append(HLA_fam)
    print(f"Sequence file output: {seq_outf}")
    print(f"Alignment file output: {aln_outf}")
    
    record_info(seq_outf)

    return 0

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="convert IMGT format to fasta format and remove redundant records")
    parser.add_argument("input", help="input IMGT alignment file")
    parser.add_argument("output_prefix", nargs="?", default=None, help="output file name prefix")
    parser.add_argument("-t", "--cut_IndexofTail", default=0, type=int, help="1-based start index of tail, which is not included in output. Trimming is based on aligned sequence, gap may be included")
    args = parser.parse_args()
    main(args.input, args.output_prefix, args.cut_IndexofTail)