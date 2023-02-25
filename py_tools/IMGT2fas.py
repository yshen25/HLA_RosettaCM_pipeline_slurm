#!/usr/bin/env python3
#converts IMGT alignment format to fasta format
#Shawn 06302020
from Bio import AlignIO, SeqIO, Seq
from pathlib import Path
import tempfile, re

def convert(infile, outfile = None):
    #--------------step 1: convert IMGT multiple sequence alignment to clustal W format-----------------
    if outfile is None:
        outfile = infile.split('.')[0] + '.fasta.nr'
    #change file name if clustal intermediate is needed
    #clustal_temp = tempfile.NamedTemporaryFile()
    clustal_temp = "test.out"
    loop_index = -1
    with open(clustal_temp, "w") as outf, open(infile, "r") as inf:
        outf.write('CLUSTAL W\n')
        #read the leading sequence as reference to decypher following sequence
        data = inf.read()
        section_pattern = re.compile(r"(?:Prot.*\n.*\n)\s\w+\*\S{2,15}\s*(.*)", re.MULTILINE)
        template = []
        for item in re.finditer(section_pattern, data):
            template.append(item.group(1).replace(' ', '').replace('.', '-'))
        inf.seek(0)
        #get the length of signal peptide
        for i in range(8): #retrieve the 9th line
            inf.readline()
        mp_start = inf.readline().find("|", 21)
        signal_peptide = inf.readline()[20:mp_start]
        sp_length = len(signal_peptide.replace(" ", ""))
        inf.seek(0)
        #convert IMGT file to clustal w
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
                    continue #skip alleles that not express
                else:
                    old_rec = line1[1].replace(' ', '')
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
    #------------------------------step 2: convert clustal file to fasta format------------------------------
    with open(clustal_temp, "r") as a, open(outfile, "w") as b:
        alignments = AlignIO.read(a, 'clustal')
        HLA_fam_last = ''
        #template_list = read_template(template_db)
        for record in alignments:
            HLA_fam_now = (''.join(record.id.replace('*', '').split(':')[0:2]))
            if not HLA_fam_now == HLA_fam_last: #and not HLA_fam_now in template_list:
                record.id = HLA_fam_now
                record.seq = Seq.Seq(str(record.seq[sp_length:])).ungap('-') # cut the signal peptide here
                SeqIO.write(record, b, 'fasta')
                HLA_fam_last = HLA_fam_now
    return outfile

#abbreviated allele name is unique id used in the pipeline
'''
def read_template(template_db):
    root_dir = Path(__file__).parent.parent.resolve()
    #template_path_abs = str(root_dir) + template_db
    #print(template_path_abs)
    templates = [record.id for record in SeqIO.parse(template_db, "fasta") ]
    return templates
    '''
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="convert IMGT format to fasta format and remove redundant records")
    parser.add_argument("input", help="input IMGT alignment file")
    parser.add_argument("output", nargs="?", default=None, help="output file name")
    args = parser.parse_args()
    print(convert(args.input, args.output))