import numpy as np

# ---------------------------------------------------------
# Constants
# ---------------------------------------------------------
SEQ_LENGTH = 1270          # Spike sequence length
NUM_SEQUENCES = 64075      # Number of sequences in BA1.fasta

# Histogram sizes
NUM_I_BINS = 55            
NUM_J_BINS = 40            

# Basic amino acids (for mutation classification)
BASIC_AA = {"K", "R", "H"}

# BA.1-specific mutation sites (0-based index) and reverse amino acids
BA1_SITES = [
    (66,  "A"),  # index 0
    (92,  "T"),  # index 1
    (139, "Y"),  # index 2
    (205, "L"),  # index 3
    (335, "G"),  # index 4
    (367, "S"),  # index 5
    (369, "S"),  # index 6
    (371, "S"),  # index 7
    (413, "K"),  # index 8
    (436, "N"),  # index 9
    (442, "G"),  # index 10
    (473, "S"),  # index 11
    (474, "T"),  # index 12
    (480, "E"),  # index 13
    (489, "Q"),  # index 14
    (492, "G"),  # index 15
    (494, "Q"),  # index 16
    (497, "N"),  # index 17
    (501, "Y"),  # index 18
    (543, "T"),  # index 19
    (610, "D"),  # index 20 (also counted in dcount)
    (651, "H"),  # index 21
    (675, "N"),  # index 22
    (677, "P"),  # index 23
    (760, "N"),  # index 24
    (792, "D"),  # index 25
    (852, "N"),  # index 26
    (950, "Q"),  # index 27
    (965, "N"),  # index 28
    (977, "L"),  # index 29
]

BA1_SITE_COUNT = len(BA1_SITES)

# ---------------------------------------------------------
# Input
# ---------------------------------------------------------
with open("BA1.fasta", "r") as f:
    s = f.read()

with open("BA1Cons.txt", "r") as ff:
    scons = ff.read()

# ---------------------------------------------------------
# Counters and histograms
# ---------------------------------------------------------
t = 0
num = 0

# Histograms of mutation counts per sequence
iicount = [0] * NUM_I_BINS           # distribution of total mutation counts
jjcount = [0] * NUM_J_BINS           # distribution of reverse mutation counts
kkcount = [[0] * NUM_J_BINS for _ in range(NUM_I_BINS)]  # joint distribution

kcount = 1
ss = []

# D614G reverse mutation count
dcount = 0

# Delta-specific mutation count
dlcount = 0

# Mutation type counts (basic vs. others)
bcount = 0
acount = 0

# Total mutation counts (over all positions and sequences)
ticount = 0
tjcount = 0

# Number of valid sequences used in analysis
scount = 0

# Mutations at BA.1 mutation sites (total vs. reverse)
mcount = np.zeros(BA1_SITE_COUNT)   # total reverse (toward ancestral) at BA.1 sites
ncount = np.zeros(BA1_SITE_COUNT)   # total mutations at BA.1 sites

# Per-position X count and mutation count
xcount = np.zeros(SEQ_LENGTH)       # number of sequences with X at each position
ycount = np.zeros(SEQ_LENGTH)       # number of sequences with a mutation at each position

# Consensus sequence as list of characters
cc = []
for i in range(SEQ_LENGTH):
    cc.append(scons[i])

# ---------------------------------------------------------
# Main loop over sequences
# ---------------------------------------------------------
for ii in range(NUM_SEQUENCES):
    # Skip header until ']'
    mark = True
    while mark:
        if s[t] == "]":
            mark = False
        t += 1

    # Read sequence characters until '>' (next header)
    mark = True
    xmark = 0
    count = 0
    bb = []

    while mark:
        if s[t] == "X":
            bb.append(s[t])
            xmark = 1
            count += 1
        elif s[t] == ">":
            mark = False
        elif s[t] != "\n":
            bb.append(s[t])
            count += 1
        t += 1

    # Count Delta-specific mutations when full length and without X
    if ii >= 0:  
        if (count == SEQ_LENGTH) and (xmark == 0):
            dlcheck = 0
            if bb[18] == "R":
                dlcheck = 1
                dlcount += 1
            if bb[52] == "G":
                dlcheck = 1
                dlcount += 1
            if bb[218] == "V":
                dlcheck = 1
                dlcount += 1
            if bb[48] == "R":
                dlcheck = 1
                dlcount += 1
            if bb[677] == "R":
                dlcheck = 1
                dlcount += 1
            if bb[946] == "N":
                dlcheck = 1
                dlcount += 1

    # Select sequences for analysis
    if count == SEQ_LENGTH:
    #if (count==SEQ_LENGTH) & (xmark==0):                     # exclude sequences with X
    #if (count==SEQ_LENGTH) & (xmark==0) & (dlcheck==0):      # exclude X or Delta-specific mutations
        scount += 1
        check = 0
        ckcount = 0

        # Replace X by consensus and count total mutations per sequence
        for j in range(SEQ_LENGTH):
            if bb[j] == "X":
                bb[j] = cc[j]
                xcount[j] += 1
            if bb[j] != cc[j]:
                ckcount += 1

        # Exclude sequences with too many mutations
        if ckcount > 25:
            check = 1
            scount -= 1

        icount = 0  # number of mutations in this sequence
        jcount = 0  # number of reverse mutations at BA.1 sites
        dd = []     # list of mutation position / amino acid (unused but kept for compatibility)
        aacount = 0 # non-basic amino acid mutation count
        bbcount = 0 # basic amino acid mutation count

        if check == 0:
            ss.append(bb)
            kcount += 1

            # Scan mutations over all positions
            for i in range(SEQ_LENGTH):
                if bb[i] != cc[i]:
                    ycount[i] += 1
                    dd.append(str(i + 1))
                    dd.append(bb[i])
                    icount += 1

                    # Classification by amino acid type
                    if bb[i] in BASIC_AA:
                        bbcount += 1
                    else:
                        aacount += 1

                    # -----------------------------
                    # Count reverse mutations at BA.1 mutation sites
                    # -----------------------------
                    for site_index, (pos, rev_aa) in enumerate(BA1_SITES):
                        if i == pos:
                            ncount[site_index] += 1  # total mutations at this BA.1 site
                            if bb[i] == rev_aa:
                                jcount += 1
                                mcount[site_index] += 1
                                # Special case: D614G reverse (position 610 in 0-based index)
                                if pos == 610:
                                    dcount += 1
                            break  # no need to check other sites for this position

        # Update histograms of mutation / reverse mutation counts
        if icount < NUM_I_BINS:
            iicount[icount] += 1
        if jcount < NUM_J_BINS:
            jjcount[jcount] += 1
        if (icount < NUM_I_BINS) and (jcount < NUM_J_BINS):
            kkcount[icount][jcount] += 1

        # Global counters
        bcount += bbcount
        acount += aacount
        ticount += icount
        tjcount += jcount

# ---------------------------------------------------------
# Output
# ---------------------------------------------------------
print("scount", scount)

print("read counts")
for i in range(SEQ_LENGTH):
    readcount = scount - xcount[i]
    print(i + 1, ",", readcount)

print("mutation counts")
for i in range(SEQ_LENGTH):
    print(i + 1, ",", ycount[i])

# Optional Check
# print("dlcount", dlcount)
# print("614D", dcount)
# print(acount, bcount)
# print(ticount, tjcount)

print("total mutations at BA.1 mutation sites")
print(
    mcount[0], ",", mcount[1], ",", mcount[2], ",", mcount[3], ",", mcount[4], ",",
    mcount[5], ",", mcount[6], ",", mcount[7], ",", mcount[8], ",", mcount[9], ",",
    mcount[10], ",", mcount[11], ",", mcount[12], ",", mcount[13], ",", mcount[14], ",",
    mcount[15], ",", mcount[16], ",", mcount[17], ",", mcount[18], ",", mcount[19], ",",
    mcount[20], ",", mcount[21], ",", mcount[22], ",", mcount[23], ",", mcount[24], ",",
    mcount[25], ",", mcount[26], ",", mcount[27], ",", mcount[28], ",", mcount[29]
)

print("reverse mutations at BA.1 mutation sites")
print(
    ncount[0], ",", ncount[1], ",", ncount[2], ",", ncount[3], ",", ncount[4], ",",
    ncount[5], ",", ncount[6], ",", ncount[7], ",", ncount[8], ",", ncount[9], ",",
    ncount[10], ",", ncount[11], ",", ncount[12], ",", ncount[13], ",", ncount[14], ",",
    ncount[15], ",", ncount[16], ",", ncount[17], ",", ncount[18], ",", ncount[19], ",",
    ncount[20], ",", ncount[21], ",", ncount[22], ",", ncount[23], ",", ncount[24], ",",
    ncount[25], ",", ncount[26], ",", ncount[27], ",", ncount[28], ",", ncount[29]
)