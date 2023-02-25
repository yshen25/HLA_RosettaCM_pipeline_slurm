#!usr/bin/env python3
#parse pssm file from psiblast, then output mtx file and checkpoint file
#Shawn 01262020

epsilon = 0.000001
ncbicodes = "*A*CDEFGHIKLMNPQRSTVWXY*****"
aafreq = [0, 0.07805, 0, 0.01925, 0.05364, 0.06295, 0.03856, 0.07377, 0.02199, 0.05142, 0.05744, 0.09019,
    0.02243, 0.04487, 0.05203, 0.04264,	0.05129, 0.07120, 0.05841, 0.06441, 0.01330, 0, 0.03216, 0, 0, 0, 0, 0]
    # original version is 26 values, two 0s are appended to the list.

import numpy as np
import re
import math

class chkparser():
    def __init__(self):
        self.len = 0
        self.seq = ""
        self.mtx_data = []
        self.chk_data = []

    def error(self, message):
        raise ValueError(message)
    
    def MatIndex(self, aa:str):
        key = {"A":0,"C":4,"D":3,"E":6,"F":13,"G":7,"H":8,"I":9,"K":11,"L":10,"M":12,"N":2,"P":14,"Q":5,"R":1,"S":15,"T":16,"V":19,"W":17,"Y":18,"X":20}
        return key[aa]

    def BLOSUM62(self, index):
        #21 aa version
        #ARNDCQEGHILKMFPSTWYVX
        matrix = np.array([[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,0],
                [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1],
                [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,-1],
                [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,-1],
                [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-2],
                [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,-1],
                [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,-1],
                [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1],
                [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,-1],
                [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-1],
                [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-1],
                [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,-1],
                [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-1],
                [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-1],
                [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2],
                [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0],
                [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,0],
                [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-2],
                [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-1],
                [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-1],
                [0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1]])
        return matrix[index]

    def parsepssm(self, pssm_file):
        with open(pssm_file, "r") as pssmhandle:
            HasTitle = 0
            HasLength = 0
            HasSeq = 0
            copy_seq = False
            raw_seq = ""

            for line in pssmhandle:# loop 1, search for titile, length, and sequence
                if not HasTitle:
                    if "PssmWithParameters" in line:
                        HasTitle = 1
                    else:
                        continue

                if not HasLength:
                    if "length" in line:
                        full_length = int(re.search(r"length (\d+),", line).group(1))
                        HasLength = 1
                    else:
                        continue

                if not HasSeq:
                    if "seq-data iupacaa" in line:
                        copy_seq = True

                if copy_seq:
                    raw_seq += line.strip()
                    if "\"\n" in line:
                        copy_seq = False
                        full_seq = re.search(r"\"(\w+)\"", raw_seq).group(1)

                        if len(full_seq) != full_length:
                            self.error("!!! sequence length doesn't match !!!")
                        HasSeq = True
                    continue
                
                if HasTitle and HasLength and HasSeq:
                    break
            
            else: #no break
                self.error("!!! file format not recognized !!!")

            copy_freq = False
            raw_freq = []

            for line in pssmhandle:# loop 2, capture freqratio
                #print(line)
                if "freqRatios" in line:# fraqratio starts
                    copy_freq = True
                    continue

                if "      }" in line:# fraqratio ends
                    copy_freq = False
                    break

                if copy_freq:
                    val, base, power = re.search(r"{ (\d+), (\d+), (-{0,1}\d+) }", line).group(1,2,3)
                    raw_freq.append(int(val)*int(base)**int(power))

            else: #no break
                self.error("!!! fraqRatio not readable !!!")

            copy_score = False
            raw_score = []

            for line in pssmhandle:# loop 3, capture score
                if "scores" in line:# score starts
                    copy_score = True
                    continue

                if "      }" in line:# score ends
                    copy_score = False
                    break

                if copy_score:
                    raw_score.append(int(line.replace(",", "").strip()))

            else: #no break
                self.error("!!! score not readable !!!")

            full_freq = np.array(raw_freq).reshape((full_length,28))
            full_score = np.array(raw_score).reshape((full_length,28))

        return full_length, full_seq, full_freq, full_score

    def GenerateContent(self, pssm_file, coverage = None):
        
        full_length, full_seq, full_freq, full_score = self.parsepssm(pssm_file)
        
        if coverage is None:
            start = 0
            stop = full_length
        else:
            start, stop = map(int, coverage.split("-"))
            start = start - 1 # input is 1-based, change to 0-based
            # stop is included in input, so no -1 needed

        self.seq = full_seq[start: stop]
        self.len = stop - start

        sxx = 0
        sxy = 0

        for i in range(0, full_length):
            for j in range(0, 28):
                if full_freq[i][j] > epsilon and aafreq[j] > epsilon:
                    
                    x = math.log(full_freq[i][j] / aafreq[j])
                    y = full_score[i][j]

                    sxx += (y * y) * x * x # Weight by y^2
                    sxy += (y * y) * x * y
        # Estimate original scaling factor by weighted least squares regression
        scale = 100.0 * sxy / sxx

        for i in range(start, stop):
            mtx_line = []
            chk_line = []
            chk_line.append(full_seq[i])

            for j in range(0, 28):
                if ncbicodes[j] == "*":
                    
                    mtx_line.append("-32768")

                else:

                    if ncbicodes[j] != "X":
                        
                        chk_line.append(f"{full_freq[i][j]:.6f}")
                        
                    if full_freq[i][j] > epsilon:

                        mtx_line.append(str(round(scale * math.log(full_freq[i][j] / aafreq[j]))))

                    else:

                        mtx_line.append(str(100*self.BLOSUM62((self.MatIndex(full_seq[i]), self.MatIndex(ncbicodes[j])))))

            self.mtx_data.append("  ".join(mtx_line) + "\n")
            self.chk_data.append(" ".join(chk_line) + "\n")        

        return 0

    def Output(self, mtx_file, chk_file):
        
        with open(mtx_file, "w") as mtx_handle:
            mtx_handle.write(f"{self.len}\n")
            mtx_handle.write(f"{self.seq}\n")
            mtx_handle.write("0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n") #12 0s
            mtx_handle.writelines(self.mtx_data)

        with open(chk_file, "w") as chk_handle:
            chk_handle.write(f"{self.len}\n")
            chk_handle.writelines(self.chk_data)
            chk_handle.write("END")

def concatenate(parser1, parser2):
    parser3 = chkparser()
    parser3.len = parser1.len + parser2.len
    parser3.seq = parser1.seq + parser2.seq
    parser3.mtx_data = parser1.mtx_data + parser2.mtx_data
    parser3.chk_data = parser1.chk_data + parser2.chk_data

    return parser3


def single(pssm_file, mtx_file, chk_file, coverage=None):
    parser = chkparser()
    parser.GenerateContent(pssm_file, coverage)
    parser.Output(mtx_file, chk_file)

    return 0

def double(pssm_file1, pssm_file2, mtx_file, chk_file, coverage1=None, coverage2=None):
    parser1 = chkparser()
    parser1.GenerateContent(pssm_file1, coverage1)
    parser2 = chkparser()
    parser2.GenerateContent(pssm_file2, coverage2)

    parser3 = concatenate(parser1, parser2)
    parser3.Output(mtx_file, chk_file)

    return 0

#single("/home/shawn/work_bench/DOWN_SYNC/test.pssm", "/home/shawn/work_bench/DOWN_SYNC/test.mtx", "/home/shawn/work_bench/DOWN_SYNC/test.chk")