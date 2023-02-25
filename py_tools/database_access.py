#!/usr/bin/env python3
#Shawn 07292020
#manipulate .csv files
#job status code: - not performed; + performed; # skipped; * low priority; ! incomplete; ? allele not found (check_status only)
import csv
import os
import shutil

def read_csv(infile):
    #read csv files into list of dict
    with open(infile, "r", newline='') as inf:
        csv_reader = csv.DictReader(inf)
        out_list = [row for row in csv_reader]
        header = csv_reader.fieldnames
    
    return out_list, header

def write_csv(out_list, header, outfile):
    # write list of dict into csv
    with open(outfile, "w", newline='') as ouf:
        csv_writer = csv.DictWriter(ouf, header)
        csv_writer.writeheader()
        csv_writer.writerows(out_list)
    
    return 0

def allele_lookup(record_list, Allele):
    #lookup by allele name in list given bu read_csv
    for row in record_list:
        if row['Allele'] == Allele:
            return row

    return

def gene_lookup(record_list, query_A, query_B=None):
    #lookup by gene name in list given by read_csv
    if query_B is None:
        for row in record_list:
            if query_A == row['A_gene']:
                return row

    else:
        for row in record_list:
            if row['A_gene'] == query_A and row['B_gene'] == query_B:
                return row
    
    return


class job_init():

    def __init__(self, job_name, target_D, record_file1, record_file2=None, filter_by_length=1, len_threshold=0):
        if record_file2 is None:
            self.MHC_class = 1
        else:
            self.MHC_class = 2
        self.job_name = job_name
        self.work_dir = target_D
        self.A_rec_file = record_file1
        self.B_rec_file = record_file2
        self.filter_by_length = filter_by_length
        self.len_threshold = len_threshold

    def len_fail(self, pep_len):
        if int(pep_len) < self.len_threshold:
            return 1
        else:
            return 0

    def arrange_target_alleles(self):
        #Class I only, generate the project file
        rec_list, _ = read_csv(self.A_rec_file)
        #Allele_name|Peptide_length|Sequence_hash
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)

        with open(f"{self.work_dir}/{self.job_name}.csv", "w", newline='') as ouf:
            pj_writer = csv.writer(ouf, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            pj_writer.writerow(["Stat", "Allele", "A_index", "A_length"])
            
            A_index = 1 #1-based, used by psiblast

            for allele in rec_list:
                A_gene = allele["Allele_name"]
                A_length = allele["Peptide_length"]
                Stat = "*"
                
                if self.filter_by_length and self.len_fail(int(A_length)):
                    Stat = "#"

                pj_writer.writerow([Stat,A_gene,A_index,A_length])

                A_index +=1

        return 0

    def combine_target_alleles(self):
        #Class II only, generate the project file  
        A_rec_list, _ = read_csv(self.A_rec_file)
        #Allele_name|Peptide_length|Sequence_hash
        B_rec_list, _ = read_csv(self.B_rec_file)

        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)

        with open(f"{self.work_dir}/{self.job_name}.csv", "w", newline='') as ouf:
            cb_writer = csv.writer(ouf, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            cb_writer.writerow(["Stat", "Allele", "A_gene", "A_index", "A_length", "B_gene", "B_index", "B_length"])
            
            A_index = 1 #1-based, used by psiblast
            B_index = 1 #1-based, used by psiblast

            for a in A_rec_list:
                A_gene = a["Allele_name"]
                A_length = a["Peptide_length"]
                for b in B_rec_list:
                    B_gene = b["Allele_name"]
                    B_length = b["Peptide_length"]

                    Stat = "*"
                
                    if self.filter_by_length and (self.len_fail(int(A_length)) or self.len_fail(int(B_length))):
                        Stat = "#"

                    cb_writer.writerow([Stat,f"{A_gene}_{B_gene}",A_gene,A_index,A_length,B_gene,B_index,B_length])

                    B_index += 1
                
                B_index = 1
                A_index +=1

        return 0

class StatusEditor():

    def __init__(self, job_file, work_dir=None):
        self.job_file = job_file
        self.result = []
        
        if work_dir is None:
            self.work_dir = os.path.dirname(job_file)
        else:
            self.work_dir = work_dir
        
        self.record_list, self.header = read_csv(job_file)
        
        if "B_gene" in self.header:
            self.MHC_class = 2
        else:
            self.MHC_class = 1

    def change_stat(self, in_row, change_to):
        in_row['Stat'] = change_to
        return in_row

    def stdout_print(self, AlleleName, Case):
        if Case == "Finish":
            Message = f"Finished: {AlleleName}"

        elif Case == "Incomplete":
            Message = f"Incomplete job: {AlleleName}"

        elif Case == "NotFound":
            Message = f"Allele not found: {AlleleName}, check query allele(s)"

        print(Message)
        
        return

    def output_exist(self, sub_dir):
        # 1=exist, 0=not
        # print(f"{self.work_dir}/{sub_dir}/models/output")
        if os.path.exists(f"{self.work_dir}/{sub_dir}/models/output"):
            return 1
        else:
            return 0

    def update_job_file(self):
        write_csv(self.record_list, self.header, f"{self.job_file}.temp")
        os.remove(self.job_file)
        shutil.move(f"{self.job_file}.temp", self.job_file)

    def through_job_file(self, status=None):
        # loop through job file and return one row or record each time
        # if status is provided, return rows that in corresponding status
        # read-only
        for row in self.record_list:
            if (row['Stat'] in status) or (status is None):
                yield row

            else:
                continue

        return

    def check_stat(self, query_list_file=None, return_list=0):
        result_list = [] # record check result

        if query_list_file is None:
            # check status for every finished record
            if self.MHC_class == 1:
                for row in self.record_list:
                    if row['Stat'] in ("+", "!"):
                        sub_dir = row['Allele']

                        if self.output_exist(sub_dir):
                            row = self.change_stat(row, "+")
                            result_list.append(row) # in check-all mode, result list only record finished jobs
                        
                        else:
                            row = self.change_stat(row, "!")
                            self.stdout_print(sub_dir, "Incomplete")
                    else:
                        continue

            if self.MHC_class == 2:
                
                for row in self.record_list:
                    
                    if row['Stat'] in ("+", "!"):
                        sub_dir = f"{row['A_gene']}/{row['B_gene']}"
                        
                        if self.output_exist(sub_dir):
                            row = self.change_stat(row, "+")
                            result_list.append(row)
                        
                        else:
                            row = self.change_stat(row, "!")
                            self.stdout_print(sub_dir, "Incomplete")
                    
                    else:
                        continue

        elif query_list_file:
            # check status for every allele on a list input
            with open(query_list_file, "r") as list_file:
                query_list = [line.strip() for line in list_file]
                
            if self.MHC_class == 1:
                
                for row in self.record_list:
                    
                    if row['Allele'] in query_list:
                        
                        if row['Stat'] in ("+", "!"):
                            sub_dir = row['Allele']
                            
                            if self.output_exist(sub_dir):
                                row = self.change_stat(row, "+")

                            else:
                                row = self.change_stat(row, "!")
                                self.stdout_print(sub_dir, "Incomplete")
                        
                        else:
                            continue
                        result_list.append(row)
                        query_list.remove(row['Allele'])

            if self.MHC_class == 2:
                
                for row in self.record_list:
                    
                    if row['Allele'] in query_list:
                        
                        if row['Stat'] in ("+", "!"):
                            sub_dir = f"{row['A_gene']}/{row['B_gene']}"

                            if self.output_exist(sub_dir):
                                row = self.change_stat(row, "+")

                            else:
                                row = self.change_stat(row, "!")
                                self.stdout_print(sub_dir, "Incomplete")

                        result_list.append(row)
                        query_list.remove(row['Allele'])
            
            if len(query_list) > 0:
                print(f"Alleles not found:\n{query_list}")

            write_csv(result_list, self.header, query_list_file+".status")

        self.update_job_file()

        if return_list:
            return result_list
        
        return

    def switch_priority(self, activate=0, query_list_file=None):
        # change job status. 0: '-' -> '*' ; 1: '*' -> '-'
        # if no list file provided, change all
        if query_list_file:
            with open(query_list_file, "r") as list_file:
                query_list = [line.strip() for line in list_file]
            result_list = []

            for row in self.record_list:
                if row['Allele'] in query_list:
                    if activate and row['Stat'] == "*":
                        row = self.change_stat(row, "-")
                        print(f"Activate job: {row['Allele']}")
                    elif (not activate) and row['Stat'] == "-":
                        row = self.change_stat(row, "*")
                        print(f"Inactivate job: {row['Allele']}")
                    result_list.append(row)
                    query_list.remove(row['Allele'])
            
            if len(query_list) > 0:
                print(f"Alleles not found:\n{query_list}")

            write_csv(result_list, self.header, query_list_file+".status")

        else:
            for row in self.record_list:
                if activate and row['Stat'] == "*":
                    row = self.change_stat(row, "-")
                elif (not activate) and row['Stat'] == "-":
                    row = self.change_stat(row, "*")
            
        self.update_job_file()
        
        return 0

    def reset(self, query_list_file=None):
        # change job status: '!' or '+' -> '-'
        # if no list file provided, change all
        if query_list_file:
            with open(query_list_file, "r") as list_file:
                query_list = [line.strip() for line in list_file]
            result_list = []

            for row in self.record_list:
                if row['Allele'] in query_list:
                    if row['Stat'] == "!":
                        row = self.change_stat(row, "-")
                        print(f"Reset job: {row['Allele']}")
                    result_list.append(row)
                    query_list.remove(row['Allele'])
            
            if len(query_list) > 0:
                print(f"Alleles not found:\n{query_list}")

            write_csv(result_list, self.header, query_list_file+".status")

        else:
            for row in self.record_list:
                if row['Stat'] == "!":
                    row = self.change_stat(row, "-")
            
        self.update_job_file()

    def launch_queue(self, numofjobs, query_list_file=None):
        # api for job_launcher
        # feed records from job file to job_launcher, then update job file
        queue = []
        if query_list_file:
            with open(query_list_file, "r") as list_file:
                query_list = [line.strip() for line in list_file]

            index = 0
            for row in self.record_list:
                
                if row['Allele'] in query_list and row['Stat'] == ("-" or "*"):
                    print(f"Queue job: {row['Allele']}")
                    queue.append((index,row))
                    query_list.remove(row['Allele'])
                
                if len(queue) >= numofjobs:
                    break

                index += 1

            if len(query_list) > 0:
                print(f"Alleles not found:\n{query_list}")

        else:
            index = 0
            for row in self.record_list:
                if row['Stat'] == "-":
                    queue.append((index, row))
                    print(f"Queue job: {row['Allele']}")
                index += 1

                if len(queue) >= numofjobs:
                    break
        
        return queue

    def queue_confirm(self, queue):
        for queue_record in queue:
            index, _ = queue_record
            self.record_list[index] = self.change_stat(self.record_list[index], "+")

        self.update_job_file()

        return 0


def initial(args):
    if len(args.record_file) == 2:
        record_file1 = args.record_file[0]
        record_file2 = args.record_file[1]
    else:
        record_file1 = args.record_file[0]
        record_file2 = None

    if args.min_length == 0:
        filter_by_length = 0
    else:
        filter_by_length = 1

    operator = job_init(args.job_name, args.output_dir, record_file1, record_file2, filter_by_length=filter_by_length, len_threshold=args.min_length)
    
    if operator.MHC_class == 1:
        operator.arrange_target_alleles()
    elif operator.MHC_class == 2:
        operator.combine_target_alleles()

    return 0

def check(args):
    operator = StatusEditor(args.job_file)
    operator.check_stat(args.query_list)
    return 0

def activate(args):
    operator = StatusEditor(args.job_file)
    operator.switch_priority(activate=1, query_list_file=args.query_list)
    return 0

def inactivate(args):
    operator = StatusEditor(args.job_file)
    operator.switch_priority(activate=0, query_list_file=args.query_list)
    return 0

def reset(args):
    operator = StatusEditor(args.job_file)
    operator.reset(query_list_file=args.query_list)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=lambda args: parser.print_help())
    subparser = parser.add_subparsers(help="choose function")

    # parser for initialization command
    parser_a = subparser.add_parser("initial", help="create job file from record file (IMGT_input.py output)")
    parser_a.set_defaults(func=initial)
    parser_a.add_argument("job_name", help="job file will be named job_name.csv")
    parser_a.add_argument("record_file", nargs='+', help="'.csv' record file(s), give 2 for class II molecules")
    parser_a.add_argument("-d","--output_dir", help="output directory", default=".")
    parser_a.add_argument("-m", "--min_length", type=int, default=100, help="inactivate alleles shorter than threshold. if set to 0, no filter will be performed")

    # parser for check function
    parser_b = subparser.add_parser("check", help="check job if completed")
    parser_b.set_defaults(func=check)
    parser_b.add_argument("job_file", help="job file")
    parser_b.add_argument("query_list", nargs="?", help="a file containing alleles for checking, one allele each line")

    # parser for activation function
    parser_c = subparser.add_parser("activate", help="change job status from '*' to '-'")
    parser_c.set_defaults(func=activate)
    parser_c.add_argument("job_file", help="job file")
    parser_c.add_argument("query_list", nargs="?", help="alleles to be activated")

    # parser for inactivation function
    parser_d = subparser.add_parser("inactivate", help="change job status from '-' to '*'")
    parser_d.set_defaults(func=inactivate)
    parser_d.add_argument("job_file", help="job file")
    parser_d.add_argument("query_list", nargs="?", help="alleles to be inactivated")

    # parser for resetting function
    parser_e = subparser.add_parser("reset", help="change incomplete job status from '!' to '-'")
    parser_e.set_defaults(func=reset)
    parser_e.add_argument("job_file", help="job file")
    parser_e.add_argument("query_list", nargs="?", help="alleles for resetting")

    args = parser.parse_args()
    args.func(args)
