#takes in a list of fasta files from spades (maybe quast) and runs blast on them,
# identifying the strains they came from (MLST)

import argparse
import os
import glob
import subprocess

def getfasta(path):
    fastalist = []
    for x in path:
        blurgh = glob.glob(x+'*.fasta')
        fastalist.append(blurgh)
    return fastalist


def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)


def mistargs(fastalist, outpath, testtype, mistfolder):
    for file in fastalist:
        strain, extension = os.path.splitext(os.path.basename(file))
        missed=('/usr/local/bin/MIST/MIST.exe', '-b',
                '-j', outpath+strain+testtype+'.json',
                '-a', mistfolder+'alleles/',
                '-t', mistfolder+testtype+'.markers',
                file)
        yield missed


def runmist(fastalist, outpath, testtype, mistfolder):
    for missed in mistargs(fastalist, outpath, testtype, mistfolder):
        subprocess.call(missed)


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outpath', default='./mistout/')
    parser.add_argument('-m', '--mistfolder', default='/home/cintiq/Desktop/campylobacterjejuni/')
    parser.add_argument('-t', '--testtype', required=True, help='type of test, ex. CGF119')
    parser.add_argument('path', nargs='+')
    return parser.parse_args()


def process(path, outpath, testtype, mistfolder):
    listlist = getfasta(path)
    for fastalist in listlist:
        pathfinder(outpath)
        runmist(fastalist, outpath, testtype, mistfolder)


def main():
    args = arguments()
    process(args.path, args.outpath, args.testtype, args.mistfolder)

if __name__ == '__main__':
    main()