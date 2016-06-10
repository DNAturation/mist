#takes input from parsed JSON files from MIST and sets threshhold
#report part1 = dictionary of strains, values of genes they are missing
#report part2 = dictionary of genes, values of strains that lack those genes

import argparse
import json
import glob
import os
import string
import multiprocessing
from functools import partial


def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)

def fileget(path):
    files = glob.glob(path+'*.json')
    return files

def reader(files):
    with open(files) as f:
        data = json.load(f)
    return data


def writer(data, testtypename, files):
    '''previously wrote data to files (hence the name), now provides the strain name,
    number of genes the genome does not have, and a list of genes that did not match for the genome
    used as primary source of information for part1 of report'''
    genes = data['Results'][0]['TestResults'][testtypename]
    genesmissing = 0
    genelist = []
    for gene in genes:
        if genes[gene]['CorrectMarkerMatch'] == False:
            genelist.append(gene)
            genesmissing+=1
    return os.path.splitext(os.path.basename(files))[0][:-len(testtypename)], genesmissing, genelist



def genes(genelist, dwriter):
    '''takes in a list of all genes that have misses and a dictionary of k=strain v=list of missing genes,
    returns a dictionary of k=genes v=strains that lack it for part 2 of the report'''
    dmisslist={}
    for gene in genelist:
        dmisslist[gene]=[]
    for strain in dwriter:
        for index, genes in enumerate(dwriter[strain]):
            dmisslist[genes].append(strain)
    return dmisslist



def JSONwriter(outpath, outfile, dwriter, dmisslist):
    '''writes the report into a json file'''
    letterlist = ['']
    letterlist = letterlist + list(string.ascii_lowercase)
    for letter in letterlist:
        if os.path.isfile(os.path.join(outpath, outfile+letter+'.json')):
            pass
        else:
            with open(os.path.join(outpath, outfile+letter+'.json'), 'a') as f:
                data = {'GenomesMissingGenes': dwriter, 'GenesMissingGenomes': dmisslist}
                json.dump(data, f, indent=4, sort_keys=True)
            return


def genetotal(testtype):
    '''returns all genes in the .markers file'''
    with open(testtype, 'r') as f:
        data = f.readlines()
        genelist = []
        for x in range(1, len(data)):
            genelist.append(data[x].split()[0])
        return genelist

def mult(item, testtypename):
    '''generates a list of dictionaries that contain strains that a gene is missing from in parallel'''
    dwriter = {}
    # d = manager.dict({})
    if not os.stat(item).st_size == 0:
        data = reader(item)
        missingno = writer(data, testtypename, item)
        # d[missingno[0]]=missingno[1]
        dwriter[missingno[0]] = missingno[2]
    else:
        print('Skipping file '+item+' due to empty .json file')
    return dwriter

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outpath', default='./mistfail/')
    parser.add_argument('--outfile', default='proteinsfailed')
    parser.add_argument('-t', '--testtype', required = True)
    parser.add_argument('-T', '--testtypename', required=True)
    parser.add_argument('-c', '--cores', default=multiprocessing.cpu_count())
    parser.add_argument('path')
    return parser.parse_args()

def process(path, outpath, outfile, testtype, testtypename, cores):
    pool = multiprocessing.Pool(int(cores))
    # manager = multiprocessing.Manager()
    pathfinder(outpath)
    files = fileget(path)
    dwriter={}
    genelist = genetotal(testtype)
    fatdwriter = pool.map(partial(mult, testtypename=testtypename), files)
    for dictionary in fatdwriter:
        for key, value in dictionary.items():
            dwriter[key]=value
            # mult(item, testtypename, dwriter, d)
            # data = reader(item)
            # missingno = writer(data, testtypename, item)
            # d[missingno[0]]=missingno[1]
            # dwriter[missingno[0]] = missingno[2]
    # pool.close()
    # pool.join()
    dmisslist=genes(genelist, dwriter)
    JSONwriter(outpath, outfile, dwriter, dmisslist)

def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.testtype, args.testtypename, args.cores)


if __name__ == '__main__':
    main()