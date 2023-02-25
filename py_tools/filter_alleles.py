#!/usr/bin/env python3
#filter fasta file from IMGT/HLA
#remove redundant records, keeps the most completed ones.
#Shawn 080819
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

def pick(infile, outfile = None):
    if outfile is None:
        outfile = infile + '.nr'
    #iterate through all records, and pick unique alleles.
    current_tag = ''
    flag = 0
    templates = read_template()
    with open(outfile, "w") as outf:
        for record in SeqIO.parse(infile, "fasta"):
            (allele, length) = record.description.split(" ")[1:3]
            if not allele.endswith(('N', 'L', 'S', 'C', 'A', 'Q')):
                allele_tag = ''.join(allele.replace('*', '').split(':')[0:2])
                length = int(length)
                if not allele_tag in templates:
                    if allele_tag != current_tag:
                        if flag == 1:
                            new_record = SeqRecord(max_seq, id=current_tag, description=max_id)
                            SeqIO.write(new_record, outf, "fasta")
                        max_id = record.id
                        max_length = length
                        max_seq = record.seq
                        current_tag = allele_tag
                        flag = 1
                    else:
                        if length > max_length:
                            max_length = length
                            max_id = record.id
                            max_seq = record.seq
                else:
                    continue
        new_record = SeqRecord(max_seq, id=current_tag, description=max_id)
        SeqIO.write(new_record, outf, "fasta")
    return outfile

#abbreviated allele name is unique id used in the pipeline
def read_template(template_seq_file = "database/blastDB/templates.fasta"):
    current_path_abs = Path(__file__).parent.parent.resolve()
    template_path_abs = str(current_path_abs.joinpath(template_seq_file))
    #print(template_path_abs)
    templates = [record.id for record in SeqIO.parse(template_path_abs, "fasta") ]
    return templates

if __name__ == "__main__":
    import sys
    try:
        pick(sys.argv[1], sys.argv[2])
    except:
        pick(sys.argv[1])