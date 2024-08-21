import numpy as np
import pandas as pd
import sys

MAX = 2000 # data size used to create a consensus sequence
MLEN = 1400 # maximum length of sequence 

aminoacid = "ACDEFGHIKLMNPQRSTVWY"

f = open('KP3_1_1.fasta', 'r') # reading sequence data
s = f.read()
f.close()

t = 0

count = np.zeros(MAX) # sequence length counter

bb = [[] for _ in range(MAX)]

for ii in range (MAX):
    mark = 1
    while mark == 1:
        if(s[t] == "]"): # skipped till the end of header
            mark = 0
        t = t + 1
    mark = 1

    # reading sequence
    while mark == 1:
        if t >= len(s):          
            sys.exit() # exit when data is no longer available
        elif s[t] == ">":
            mark = 0
#            totalcount += 1
        elif s[t] != "\n":
            bb[ii] += s[t]
            count[ii] += 1
        t = t + 1

histo = np.zeros(MLEN)
for ii in range(MAX):
    for jj in range(MLEN):
        if count[ii] == jj:
            histo[jj] += 1

max = 0
for ii in range(MLEN):
    if (histo[ii]>max):
        max = histo[ii]
        SLEN = ii

seq = ""
for ss in range(SLEN):
    countaa= np.zeros(SLEN)
    for ii in range (MAX):
        if count[ii]==SLEN:
            for jj in range(20):
                if(bb[ii][ss] == aminoacid[jj]):
                    countaa[jj] += 1
    max = 0
    for ii in range(20):
        if (countaa[ii]>max):
            max = countaa[ii]
            argmax = ii
    seq += aminoacid[argmax]
print(seq)

