#reads the report file made by mistdistrept, and symlinks fasta files that were run by mist
# that have below a threshold (set by user) amount of missing genes

import json
import os
import argparse

def pathfinder(outpath):
    '''makes a directory to store symlinked files'''
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)

def reader(path):
    '''opens the report file made by mistdistrept and stores it into data variable'''
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def countergenes(data):
    '''generates a list of the strains present and the amount of genes that they are missing'''
    strainlist = data['GenomesMissingGenes']
    for strain in strainlist:
        missingno = len(strainlist[strain])
        yield missingno, strain


def filter(missingno, strain, threshhold, passlist):
    '''determines which strains pass the cutoff, and appends them to a list that is created outside
    this function within the process function'''
    if missingno <= threshhold:
        passlist.append(strain)

def actor(passlist, outpath, mistout, testtype):
    '''does the symlinking based on the list of passed strains provided by the filter function'''
    for strain in passlist:
        if strain not in os.listdir(outpath):
            os.symlink(os.path.join(mistout, strain)+'{}.json'.format(testtype), outpath+strain)

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outpath', default='./passed/')
    parser.add_argument('-thresh', '--threshhold', type=int, required=True)
    parser.add_argument('-t', '--testtype', required=True)
    parser.add_argument('-m', '--mistout', default='/home/cintiq/PycharmProjects/misty/mistout/')
    parser.add_argument('path')
    return parser.parse_args()

def process(path, outpath, threshhold, mistout, testtype):
    pathfinder(outpath)
    data = reader(path)
    genomepasslist = []
    missingno = countergenes(data)
    for misses, strain in missingno:
        filter(misses, strain, threshhold, genomepasslist)
        actor(genomepasslist, outpath, mistout, testtype)

def main():
    args = arguments()
    process(args.path, args.outpath, args.threshhold, args.mistout, args.testtype)

if __name__ == '__main__':
    main()