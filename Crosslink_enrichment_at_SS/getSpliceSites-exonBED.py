'''
Created on Aug 20, 2018

The script will return start and end positions of the input BED file.
'''

import sys

def get_ss(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    foutEnd = open(fout_fname + "-5SS.bed", "w")
    foutStart = open(fout_fname + "-3SS.bed", "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        pos1 = col[1]
        pos2 = col[2]
        cDNA = ""
        strand = col[5]
        
        if strand == "+":
            start = int(pos1)
            end = int(pos2)
        elif strand == "-":
            start = int(pos2)
            end = int(pos1)
            
        foutStart.write(chr + '\t' + str(start) + '\t' + str(start+1) + '\t' + cDNA + '\t' + "" + '\t' + strand + '\n')
        foutEnd.write(chr + '\t' + str(end) + '\t' + str(end+1) + '\t' + cDNA + '\t' + "" + '\t' + strand + '\n')
        line = fin.readline()
    foutEnd.close()
    foutStart.close()

if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    get_ss(fin_fname, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python2 getSpliceSites-exonBED.py input_fname.bed output_fname"
