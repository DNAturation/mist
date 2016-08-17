# runs the MIST main program. takes in a list of fasta files from spades (maybe quast) and runs blast on them,
# identifying the strains they came from (MLST)

import argparse
import os
import glob
import subprocess
import multiprocessing
import shutil
import csv
import json
import re

def getfasta(path):
    '''takes in multiple directories, returns a list of all fasta files within a directory, in a bigger list'''
    fastalist = []
    for direc in path:  # path is a list due to nargs+ argument
        files = glob.glob(os.path.join(direc, '*.fasta'))
        for file in files:
            fastalist.append(file)
    return fastalist


def pathfinder(outpath):
    '''makes output directory as well as a temporary directory for the temp files MIST generates'''
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)
    if not os.access(outpath+'temp/', os.F_OK):
        os.mkdir(outpath+'temp/')

def testnamegetter(testtype):
    with open(testtype, 'r') as f:
        try: #try accessing markers file as .json file first
            data = json.load(f)
            for genome, keys in data.items():
                for key in keys:
                    if re.match('T(est)?\.?[-\._ ]?Name.*', key, flags=re.IGNORECASE):
                        return keys[key]

        except ValueError: #if access as .json file fails, try to access as csv file
            f.seek(0)
            reader=csv.reader(f, delimiter='\t')
            next(reader, None)
            for x in reader:
                testname=x[1]
                return testname



def mistargs(fastalist, outpath, testtypename, testtype, alleles):
    '''all arguments necessary for running MIST'''
    for file in fastalist:
        strain, extension = os.path.splitext(os.path.basename(file))
        if not os.path.isfile(outpath+strain+testtypename+'.json'): #skips file if the json for it already exists(has already been run on)
            if not os.stat(file).st_size == 0: #skips empty fasta files
                # missed=('mist', '-b',
                missed = ('mono', '/home/phac/kye/assemblies_for_ed/Release/MIST.exe', '-b',
                        '-j', outpath+strain+testtypename+'.json',
                        '-a', alleles,
                        '-t', testtype,
                        '-T', outpath+'temp/'+strain+'/',
                        file)
                yield missed, strain
            else:
                print('skipping strain '+strain+' due to .fasta file being an empty file')
        else:
            print('skipping strain '+strain+' due to .json file for this test already existing')


def runmist(missed, outpath, strain):
    '''calls MIST for each of the arguments provided'''
    if not os.access(outpath+'temp/'+strain+'/', os.F_OK): #makes a specific folder for each strain for their temp files
        os.mkdir(outpath+'temp/'+strain+'/')
    subprocess.call(missed)
    shutil.rmtree(outpath+'temp/'+strain+'/') #removes the strain specific temp folder and everything inside


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outpath', default='./mistout/')
    parser.add_argument('-a', '--alleles', default='alleles/')
    parser.add_argument('-t', '--testtype', required=True, help='path to and type of test/markers file, ex. CGF119')
    parser.add_argument('-c', '--cores', default=multiprocessing.cpu_count(), help='number of cores to run on')
    parser.add_argument('path', nargs='+')
    return parser.parse_args()


def process(path, outpath, testtype, alleles, cores):
    fastalist = getfasta(path)
    pool = multiprocessing.Pool(int(cores))
    testtypename=testnamegetter(testtype)
    # for fastalist in listlist:#goes through each individual list of fastas within the list of lists
    pathfinder(outpath)
    margs = mistargs(fastalist, outpath, testtypename, testtype, alleles)
    for missed, strain in margs: #run MIST in parallel
        pool.apply_async(runmist, args=(missed, outpath, strain))
    pool.close()
    pool.join()


def main():
    args = arguments()
    process(args.path, args.outpath, args.testtype, args.alleles, args.cores)

if __name__ == '__main__':
    main()