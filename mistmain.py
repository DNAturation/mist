#takes in a list of fasta files from spades (maybe quast) and runs blast on them,
# identifying the strains they came from (MLST)

import argparse
import os
import glob
import subprocess
import multiprocessing
import shutil

def getfasta(path):
    '''takes in multiple directories, returns a list of all fasta files within a directory, in a bigger list'''
    fastalist = []
    for x in path:
        blurgh = glob.glob(x+'*.fasta')
        fastalist.append(blurgh)
    return fastalist


def pathfinder(outpath):
    '''makes output directory as well as a temporary directory for the temp files MIST generates'''
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)
    if not os.access(outpath+'temp/', os.F_OK):
        os.mkdir(outpath+'temp/')


def mistargs(fastalist, outpath, testtypename, testtype, alleles):
    '''all arguments necessary for running MIST'''
    for file in fastalist:
        strain, extension = os.path.splitext(os.path.basename(file))
        if not os.path.isfile(outpath+strain+testtypename+'.json'): #skips file if the json for it already exists(has already been run on)
            if not os.stat(file).st_size == 0: #skips empty fasta files
                # missed=('/home/cintiq/Desktop/campylobacterjejuni/mist/bin/Release/MIST.exe', '-b',
                missed=('mist', '-b',
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
    parser.add_argument('-T', '--testtypename', required=True, help='name of test run')
    parser.add_argument('-t', '--testtype', required=True, help='path to and type of test/markers file, ex. CGF119')
    parser.add_argument('-c', '--cores', default=multiprocessing.cpu_count(), help='number of cores to run on')
    parser.add_argument('path', nargs='+')
    return parser.parse_args()


def process(path, outpath, testtypename, testtype, alleles, cores):
    listlist = getfasta(path)
    pool = multiprocessing.Pool(int(cores))
    for fastalist in listlist:#goes through each individual list of fastas within the list of lists
        pathfinder(outpath)
        margs = mistargs(fastalist, outpath, testtypename, testtype, alleles)
        for missed, strain in margs: #run MIST in parallel
            pool.apply_async(runmist, args=(missed, outpath, strain))
    pool.close()
    pool.join()


def main():
    args = arguments()
    process(args.path, args.outpath, args.testtypename, args.testtype, args.alleles, args.cores)

if __name__ == '__main__':
    main()