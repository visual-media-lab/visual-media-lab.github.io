import numpy as np
import pandas as pd
import sys

MAX = 100000 # max data size
SLEN = 1270 # length of consensus sequence 
NMUT = 31 # number of mutations

f = open('BA1_1_18.fasta', 'r') # reading sequence data
s = f.read()
f.close()

fc = open('BA11Cons.txt', 'r') # reading consensus sequence
scons = fc.read()
fc.close()

lst = pd.read_csv("BA11Mut.csv").values.tolist() # reading mutation list

print("Accession No.", end = "")
for i in range (NMUT):
    print(",", lst[i][0], end = "")
print("")

t = 0
iicount = np.zeros(NMUT) # mutation count
jjcount = np.zeros(NMUT) # reverse mutation count
ss = list()
an = np.zeros(8) # accessoin number
bn = s[1:9] # accession number (tmporary)
totalcount = 1 # total count of sequences
# xcount = np.zeros(SLEN) 
# mcount = np.zeros(SLEN)
xxcount = 0
xxcount2 = 0
ccmarkcount = 0
ccunmarkcount = 0

cc = [] # consensus sequence
for i in range(SLEN):
    cc += scons[i] # consensus sequence stored in cc

for ii in range (MAX):
    mark = 1
    while mark == 1:
        if(s[t] == "]"): # skipped till the end of header
            mark = 0
        t = t + 1
    mark = 1
    k = 0
    xmark = 0 # marked 1 if X is included in the sequence
    count = 0 # sequence length counter
    bb = []
    an = bn #Accession number preserved

    # reading sequence and next header information
    while mark == 1:
        if t >= len(s):          
#            print("counts of total data", totalcount)
#            print("counts of sequences including X (different length allowed)", xxcount)
#            print("counts of sequences including X", xxcount2)
#            print("counts of sequences with X at mutational points", ccmarkcount)
#            print("counts of data for analysis", ccunmarkcount)
            sys.exit() # exit when data is no longer available
        elif s[t] == "X":
            bb += s[t]
            xmark = 1
            count += 1
        elif s[t] == ">":
            mark = 0
            bn = s[t+1:t+9]
            totalcount += 1
        elif s[t] != "\n":
            bb += s[t]
            count += 1
        t = t + 1

  # sequence matching      
    if(xmark == 1):
        xxcount += 1 # counting up the number of data including X

    if(count == SLEN):
        if(xmark == 1):
            xxcount2 += 1 # counting up the number of data including X and length is the same
        ccmark = 0
        for kk in range(NMUT):
            if(bb[lst[kk][2]-1] == "X"): # X in mutation points
                ccmark = 1
        if ccmark == 1:
            ccmarkcount += 1
        else:
            ccunmarkcount += 1
        for j in range(SLEN):
            if bb[j] == "X":
                bb[j] = cc[j] # X is replaced by consensus amino acid
#                xcount[j] += 1 # counting X in each amino acid
#            if bb[j] != cc[j]:
#                mcount[j] += 1 # counting mutation in each amino acid
        icount = 0
        jcount = 0
        mt = np.zeros(NMUT)
        if ccmark == 0: # if no XX in mutation points
            ss.append(bb)
            for i in range(SLEN):
                if bb[i] != cc[i]: # if the i-th amino acid is different from the consensus amino acid
                    icount += 1 # counting up the number of mutations

                    # checking reversions in the following sequence
                    for kk in range (NMUT):
                        if(i+1 == lst[kk][2])&(bb[i] == lst[kk][1]):
                            jcount += 1 # counting up the number of reverse mutations
                            mt[kk] = 1 # reverse mutation marker

        # print pure reversion mutant pattern
        if(icount > 0) & (icount == jcount): #if mutations are included and all mutations are reversions
            print(an, end = "")
            for kk in range (NMUT):
                print(",", mt[kk], end = "")
            print("")