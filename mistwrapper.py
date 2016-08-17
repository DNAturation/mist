#main script that runs mist and generates a report, as well as symlinks based on cutoffs

import mistprog
import mistfastasym
import mistdistrept
import update_definitions
import mistgenefilter
import mistreptfilter
import mistalleledistance
import argparse
import os
import glob
import multiprocessing
import csv
import json
import re
import time


def mistP(path, outpath, testtype, alleles, cores):
    '''runs mist'''
    mistprog.process(path, outpath, testtype, alleles, cores)

def mistF(path, fastaoutpath, threshhold, outpath, testtypename, cores):
    '''runs a userdefined threshold cutoff and symlinks fasta files that pass'''
    mistfastasym.process(path, fastaoutpath, threshhold, outpath, testtypename, cores)

def mistD(path, outpath, outfile, testtype, cores):
    '''generates a report'''
    mistdistrept.process(path, outpath, outfile, testtype, cores)

def mistU(alleles, jsons, testtypename):
    '''runs updater on mist main output, dillon's script'''
    update_definitions.process(alleles, jsons, testtypename)

def mistG(path, threshhold, testtype, markerout):
    '''makes a temporary .markers file based on genes that fall under a threshold of not being present in genomes'''
    mistgenefilter.process(path, threshhold, testtype, markerout)

def mistR(path, outpath, outfile, marker, testtypename):
    '''using the new markers file, removes genes that are no longer present from the report file'''
    mistreptfilter.process(path, outpath, outfile, marker, testtypename)

def mistA(path, outpath, outfile, stout, thresh):
    '''generates a small report on how different the strains are based on which alleles were matched in MIST'''
    mistalleledistance.process(path, outpath, outfile, stout, thresh)

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


def arguments():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subfunction')
    a_parser = subparsers.add_parser('Generate')
    a_parser.add_argument('-o', '--outpath', default='/home/cintiq/PycharmProjects/misty/mistout/', help='output folder for json files from mist')
    a_parser.add_argument('-t', '--testtype', required=True, help='path to test markers file')
    a_parser.add_argument('--distoutpath', default='./mistreport/', help='output folder for the report summary')
    a_parser.add_argument('--distoutfile', default='report', help='name of the output report file')
    a_parser.add_argument('-a', '--alleles', required=True, help='directory that contains all allele files referenced by markers file')
    a_parser.add_argument('-c', '--cores', default=multiprocessing.cpu_count(), help='number of cores to run on')
    a_parser.add_argument('--seqtyp', default='ST.csv', help='filename for the output sequence-typing report')
    a_parser.add_argument('--distanceout', default='alleledistance.csv', help='filename for the output allele distance report')
    a_parser.add_argument('--thresh', default=0, type=int, help='how many allele call differences are tolerated before assigning to different strains')

    a_parser.add_argument('path', nargs='+', help='directories of fasta files to run mist on')

    b_parser = subparsers.add_parser('Refine')
    b_parser.add_argument('--genethreshhold', type=int, required=True, help='maximum amount of missed genes tolerated to allow a strain to pass')
    b_parser.add_argument('--genomethreshhold', type=int, required=True, help='maximum number of genomes that a gene is not in before being cut off')
    b_parser.add_argument('-s', '--symlinkoutpath', default='./mistpass/', help='folder to place fasta file symlinks for those that pass the threshold')
    b_parser.add_argument('-t', '--testtype', required=True, help='test markers file eg. MLST.markers')
    b_parser.add_argument('--markerout', default='/home/cintiq/PycharmProjects/misty/marker/', help='directory to place new marker in')
    b_parser.add_argument('--reptoutpath', default='./mistreport/', help='output folder for the report summary')
    b_parser.add_argument('--reptoutfile', default='refinedreport.json', help='name of the report file')
    b_parser.add_argument('-o', '--outpath', default='/home/cintiq/PycharmProjects/misty/mistout/', help='output folder for json files from mist')
    b_parser.add_argument('-c', '--cores', default=multiprocessing.cpu_count(), help='number of cores to run')
    b_parser.add_argument('path', help='path to report file')
    return parser.parse_args()

def processGEN(path, outpath, testtype, alleles, distoutpath, distoutfile, distanceout, seqtyp, thresh, cores):
    '''Runs MIST and creates a report file. Also updates alleles with new found alleles'''
    print('Running MIST...')
    mistP(path, outpath, testtype, alleles, cores)
    print('Performing update...')
    jsonlist = glob.glob(outpath+'*./json')
    testtypename=testnamegetter(testtype)
    for json in jsonlist:
        mistU(alleles, json, testtypename)
    print('Making report...')
    mistD(outpath, distoutpath, distoutfile, testtype, cores)
    mistA(outpath, distoutpath, distanceout, seqtyp, thresh)
    print('Generate has completed')


def processREF(path, symlinkoutpath, genethreshhold, genomethreshhold, testtype,
               markerout, outpath, reptoutpath, outfile, cores):
    '''
    uses a user-defined threshold to decide which files pass. Files that pass are symlinked.
    Makes a new markers file based on the genes that pass the threshold, and generates a new MIST report file
    (simulating running MIST) based on the new markers file by removing all genes that are no longer in markers
    '''
    testtypename=testnamegetter(testtype)
    print('Symlinking passing fasta files...')
    mistF(path, symlinkoutpath, genethreshhold, outpath, testtypename, cores)
    print('Filtering genes...')
    mistG(path, genomethreshhold, testtype, markerout)
    print('Revising report based on new genes...')
    mistR(path, reptoutpath, outfile, markerout, testtypename)
    print('Refine has completed')


def main():
    args = arguments()
    if args.subfunction == 'Generate':
        processGEN(args.path, args.outpath, args.testtype, args.alleles,
                   args.distoutpath, args.distoutfile, args.distanceout, args.seqtyp, args.thresh, args.cores)
    if args.subfunction == 'Refine':
        processREF(args.path, args.symlinkoutpath, args.genethreshhold, args.genomethreshhold,
                   args.testtype, args.markerout, args.outpath, args.reptoutpath,
                   args.reptoutfile, args.cores)


if __name__ == '__main__':
    main()
