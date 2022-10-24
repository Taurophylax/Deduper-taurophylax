import re
import os
import argparse as ap

parser = ap.ArgumentParser(description='taurophylax_deduper.py ')
parser.add_argument('-u', '--umifile', help='UMI List', nargs=1, type=str)
parser.add_argument('-f', '--infile', help='SORTED SAM Input File (in.sam)', nargs=1, type=str)
parser.add_argument('-o', '--outfile', help='SORTED SAM Output File (out.sam)', nargs=1, type=str)
args = parser.parse_args()
umifile = open(args.umifile, "r")
samfile = open(args.infile, "r")
outfile = open(args.outfile, "w")

def data_points(line: str) -> list: #[umi, strand, chrom, startpos, cigar, seq]
    line = re.split(r'\t+', line) #split line by tabs
    umi = line[0][-8::]  #get uni from qname
    flag = int(line[1])
    strand = "+"
    if ((flag & 16) == 16): 
        strand = "-"
    chrom = line[2]
    startpos = int(line[3])
    cigar = line[5]
    seq = line[9]
    return [umi, strand, chrom, startpos, cigar, seq]

def find_start(cigar: str, startpos: int, strand: str) -> int:    #parse cigar string and find start position
    if (strand == "+"):
        S = re.findall(r'(\d+)S', cigar)    # Find all S
        if S: S = int(S[0])                 # Convert it to int and also drop any other S occurances. LEFT clip.
        if S: startpos -= S                 # Sub S from our start position 
    else: #if strand == "-"
        S = re.findall(r'(\d+)S', cigar)    # Find all S
        if S: S = int(S[1])                 # Convert it to int and also drop any other S occurances. RIGHT clip.
        if S: startpos += S                 # Add S to our start position 
        M = re.findall(r'(\d+)M', cigar)    # Find all M
        if M: M = [int(i) for i in M]       # Convert list.str to list.int
        if M: startpos += sum(M)            # Add all Ms to position
        D = re.findall(r'(\d+)D', cigar)    # Find all D
        if D: D = [int(i) for i in D]             
        if D: startpos += sum(D)  
        N = re.findall(r'(\d+)N', cigar)    # Find all D
        if N: N = [int(i) for i in N]             
        if N: startpos += sum(N)  
    return startpos

master_umi_list = []
for a in umifile: master_umi_list.append(a.strip()) #grab umi's from STL96.txt

chunk = [] #chunk of lines for parsing (FIFO)
dupes = 0 #duplicate counter
badumi = 0 #invalid UMIs

while 1: 
    line = samfile.readline()
    if line == "":  #if EOF, break loop
        samfile.close()
        print("Done")
        break
    elif line[0] == "@":  #if header line, add directly to file
        outfile.write(line)
    else:
        linedata = data_points(line) #Returns [umi, strand, chrom, startpos, cigar, seq]
        if linedata[0] in master_umi_list:
            startpos = find_start(linedata[4], linedata[3], linedata[1]) #pass (cigar, startpos, strand) returns adjusted startpos
            linedata[3] = startpos  #update linedata with adjusted start position
            check = [linedata[0], linedata[1], linedata[2], linedata[3]]
            if check not in chunk:
                chunk.append(check)
                if len(chunk) > 100:  #if chunk is over 100 entries, remove the oldest (FIFO)
                    chunk.pop(0)
                outfile.write(line)
            else:
                dupes += 1
        else:
            badumi += 1

print("Duplicates found: " + str(dupes))
print("Invalid UMIs: " + str(badumi))