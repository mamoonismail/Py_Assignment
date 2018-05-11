#!/usr/bin/env python3

import glob
import pandas
import matplotlib.pyplot as plt
import sys

f = open('nd2.fasta', 'r') 
sequence = f.readlines()
len(sequence)

to separate one fasta file into multiple fasta files with the filename = the header of the sequence
infile = open('nd2.fasta', 'r')
outfile = []

for line in infile:
    if line.startswith(">"):
        if (outfile != []): outfile.close()
        genename = line.strip().split('>')[1]
        filename = genename+".fasta"
        outfile = open(filename,'w')
        outfile.write(line)
    else:
        outfile.write(line)
outfile.close()

#For observed kmers
f = open('Adelomyia_melanogenys.fasta', 'r')
fasta = f.readlines()
for k in range(1,len(fasta[1])):
    counter = {}
    for line_num,line in enumerate(fasta):   #go through each line of the fastq file
        if line_num == 1:    #== is checking, = is assignment
            line = line.rstrip() #remove newline \n
            for i,base in enumerate(line[:-k+1]): #no longer count base individually, slice it in pairs, 
                kmer = line[i:i+k]  #remember not inclusive
                if kmer in counter:
                    counter[kmer] += 1  #if kmer has already in the dictionary, add it
                else:
                    counter[kmer] = 1   #if kmer is not in dictionary, then define it
    print(k,':',len(counter))
	
f = open('Adelomyia_melanogenys.fasta', 'r')
fasta = f.readlines()
len(fasta[1])

for k in range(1,len(fasta[1])+1):
                if 4**k < len(fasta[1]):
                    pos = 4**k
                else:
                    pos = len(fasta[1]) - k + 1
                print(k,':',pos)
				
