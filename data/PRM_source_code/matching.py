# coding: utf-8

import numpy as np
import sys

fw = open('WuhanS.txt', 'r') # reading original Wuhan sequence
X = fw.read()
fw.close()

fm = open('BA1Cons.txt', 'r') # reading consensus sequence
Y = fm.read()
fm.close()

print("Original Num, Original aa, Num after Mutation")

N = len(X)
M = len(Y)

gap_penalty = -2

H = np.empty((N+1,M+1), dtype='int16')
L = np.zeros((N+1,M+1), dtype='int8')

H[0,0]=0
for j in range(1,M+1):
    H[0,j] = H[0,j-1] + gap_penalty
    L[0,j] = 0
    
for i in range(1,N+1):
    H[i,0] = H[i-1,0] + gap_penalty
    L[i,0] = 2

# Horizontal:0 Diagonal:1 Vertical:2
s = np.array([0,0,0],dtype='int')

for i in range(1,N+1):
    for j in range (1,M+1):
        s[0] = H[i,j-1] + gap_penalty
        
        if X[i-1]==Y[j-1]:
            score = +1
        else:
            score = -1

        s[1] = H[i-1,j-1] + score
        s[2] = H[i-1,j] + gap_penalty
        
        H[i,j]=np.max(s)
        L[i,j]=np.argmax(s)

i=N
j=M
while i!=0 or j!=0:
    if L[i,j]==0:
        j=j-1
    elif L[i,j]==1:
        if X[i-1] != Y[j-1]:
            print(i, ",", X[i-1], ",", j)
        i=i-1
        j=j-1
    else:
        i=i-1