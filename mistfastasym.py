#reads the report file made by mistdistrept, and symlinks fasta files that were run by mist
# that have below a threshold (set by user) amount of missing genes

import json
import os
import argparse
import multiprocessing

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

def acter(passlist, outpath, mistout, testtypename):
    '''does the symlinking based on the list of passed strains provided by the filter function'''
    for strain in passlist:
        if strain not in os.listdir(outpath):
            os.symlink(os.path.abspath(os.path.join(mistout, strain)+'{}.json'.format(testtypename)), outpath+strain+'.json')

def mult(misses, strain, threshhold, genomepasslist, outpath, mistout, testtypename):
    filter(misses, strain, threshhold, genomepasslist)
    acter(genomepasslist, outpath, mistout, testtypename)

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outpath', default='./passed/')
    parser.add_argument('-thresh', '--threshhold', type=int, required=True)
    parser.add_argument('-t', '--testtypename', required=True)
    parser.add_argument('-m', '--mistout', default='/home/cintiq/PycharmProjects/misty/mistout/')
    parser.add_argument('-c', '--cores', default=multiprocessing.cpu_count())
    parser.add_argument('path')
    return parser.parse_args()

def process(path, outpath, threshhold, mistout, testtypename, cores):
    pathfinder(outpath)
    data = reader(path)
    genomepasslist = []
    missingno = countergenes(data)
    pool=multiprocessing.Pool(int(cores))
    for misses, strain in missingno:
        pool.apply_async(mult, args=(misses, strain, threshhold, genomepasslist, outpath, mistout, testtypename))
    pool.close()
    pool.join()

def main():
    args = arguments()
    process(args.path, args.outpath, args.threshhold, args.mistout, args.testtypename, args.cores)

if __name__ == '__main__':
    main()