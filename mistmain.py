#takes in a list of fasta files from spades (maybe quast) and runs blast on them,
# identifying the strains they came from (MLST)

import argparse
import os
import glob
import subprocess
import multiprocessing
import shutil

def getfasta(path):
    fastalist = []
    for x in path:
        blurgh = glob.glob(x+'*.fasta')
        fastalist.append(blurgh)
    return fastalist


def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)
    if not os.access(outpath+'temp/', os.F_OK):
        os.mkdir(outpath+'temp/')


def mistargs(fastalist, outpath, testtypename, testtype, alleles):
    for file in fastalist:
        strain, extension = os.path.splitext(os.path.basename(file))
        if not os.path.isfile(outpath+strain+testtypename+'.json'):
            # missed=('/home/cintiq/Desktop/campylobacterjejuni/mist/bin/Release/MIST.exe', '-b',
            missed=('mist', '-b',
                    '-j', outpath+strain+testtypename+'.json',
                    '-a', alleles,
                    '-t', testtype,
                    '-T', outpath+'temp/'+strain+'/',
                    file)
            yield missed, strain
        else:
            print('skipping strain '+strain+' due to .json file for this test already existing')


def runmist(missed, outpath, strain):
    if not os.access(outpath+'temp/'+strain+'/', os.F_OK):
        os.mkdir(outpath+'temp/'+strain+'/')
    subprocess.call(missed)
    shutil.rmtree(outpath+'temp/'+strain+'/')






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
    for fastalist in listlist:
        pathfinder(outpath)
        margs = mistargs(fastalist, outpath, testtypename, testtype, alleles)
        for missed, strain in margs:
            pool.apply_async(runmist, args=(missed, outpath, strain))
    pool.close()
    pool.join()


def main():
    args = arguments()
    process(args.path, args.outpath, args.testtypename, args.testtype, args.alleles)

if __name__ == '__main__':
    main()