import numpy as np
import pandas as pd
import sys

FASTA = "BA_2_86_1.fasta"
CONS = "BA2861Cons.txt"

MAX = 200 # data size used to create a consensus sequence
MLEN = 32000 # maximum length of sequence 

aminoacid = "ATGC"

f = open( FASTA , 'r') # reading sequence data
s = f.read()
f.close()

t = 11

count = np.zeros(MAX) # sequence length counter

bb = [[] for _ in range(MAX)]

for ii in range (MAX):
    mark = 1
    # reading sequence
    while mark == 1:
        if t >= len(s):          
            sys.exit() # exit when data is no longer available
        elif s[t] == ">":
            mark = 0
            t = t + 10
        elif s[t] != "\n":
            bb[ii] += s[t]
            count[ii] += 1
        t = t + 1

#    if ii == 0:
#        with open( "0.txt" , "w") as o:
#            print(bb[0], file=o)
#    elif ii == 1:
#        with open( "1.txt" , "w") as o:
#            print(bb[1], file=o)
#    elif ii == 2:
#        with open( "2.txt" , "w") as o:
#            print(bb[2], file=o)

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
            for jj in range(4):
                if(bb[ii][ss] == aminoacid[jj]):
                    countaa[jj] += 1
    max = 0
    for ii in range(4):
        if (countaa[ii]>=max):
            max = countaa[ii]
            argmax = ii
    seq += aminoacid[argmax]

with open( CONS , "w") as o:
    print(seq, file=o)

print(SLEN)