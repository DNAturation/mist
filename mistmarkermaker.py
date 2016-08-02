# takes in a directory of fasta files to serve as a reference for genes,
# and creates a marker file for those files for use in MIST
# Note: this assumes that the marker names are the same as the file names


import json
import argparse
import os
import glob
import re
from Bio import SeqIO


def fileget(path):
    # retrieves all fasta files from given allele directory to make a .markers pointer file for MIST
    files = glob.glob(os.path.join(path, "*.fasta"))
    return files


def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)


def basename(file):
    name = os.path.splitext(os.path.basename(file))[0]
    return name


def dmarker(d, file, testname, testtype, forward, reverse, ampsize, amprange, repeat):
    # generates the dictionary for the markers file, which will be written in json format
    name = basename(file)
    d[name] = {'Test Name': testname,
               'Test Type': testtype,
               'Forward Primer': forward,
               'Reverse Primer': reverse,
               'Amplicon Size': ampsize,
               'Amplicon Range Factor': amprange,
               'Allelic Database Filename': os.path.basename(file),
               'Repeat Size': repeat}


def jsonwriter(d, outpath, testname):
    # creates the actual .markers file and writes the generated dictionary to it from dmarker
    with open(outpath+testname+'.markers', 'w') as f:
        json.dump(d, f, indent=4, sort_keys=True)


def filechecker(file):
    # extra optional thing, checks the fasta files for non-nucleotides and gives warnings if found
    with open (file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            nonNT = re.findall('[^NACTG-]+', str(record.seq))
            for mismatch in nonNT:
                print("Warning: non-DNA nucleotide '"+mismatch+"' found in "+file)


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--testname', required=True)
    parser.add_argument('--testtype', default='1')
    parser.add_argument('--forward', default='')
    parser.add_argument('--reverse', default='')
    parser.add_argument('--ampsize', default='-1')
    parser.add_argument('--amprange', default='0')
    parser.add_argument('--repeat', default='0')
    parser.add_argument('--outpath', default='./')
    parser.add_argument('path')
    return parser.parse_args()


def process(path, outpath, testname, testtype, forward, reverse, ampsize, amprange, repeat):
    files = fileget(path)
    pathfinder(outpath)
    d = {}
    for file in files:
        filechecker(file)
        dmarker(d, file, testname, testtype, forward, reverse, ampsize, amprange, repeat)
    jsonwriter(d, outpath, testname)


def main():
    args = arguments()
    process(args.path, args.outpath, args.testname,
            args.testtype, args.forward, args.reverse,
            args.ampsize, args.amprange, args.repeat)

if __name__ == "__main__":
    main()
