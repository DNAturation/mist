#main script that runs mist and generates a report, as well as symlinks based on cutoffs
#named spray because spraying water creates a mist

import mistmain
import mistfastasym
import mistdistrept
import update_definitions
import mistgenefilter
import mistreptfilter
import argparse
import os
import glob
import multiprocessing

def mistM(path, outpath, testtypename, testtype, alleles, cores):
    '''runs mist'''
    mistmain.process(path, outpath, testtypename, testtype, alleles, cores)

def mistF(path, fastaoutpath, threshhold, outpath, testtypename):
    '''runs a userdefined threshold cutoff and symlinks fasta files that pass'''
    mistfastasym.process(path, fastaoutpath, threshhold, outpath, testtypename)

def mistD(path, outpath, outfile, testtype, testtypename):
    '''generates a report'''
    letter = mistdistrept.process(path, outpath, outfile, testtype, testtypename)
    return letter

def mistU(alleles, jsons, testtypename):
    '''runs updater on mist main output, dillon's script'''
    update_definitions.process(alleles, jsons, testtypename)

def mistG(path, threshhold, testtype, testtypename, markerout):
    '''makes a temporary .markers file based on genes that fall under a threshold of not being present in genomes'''
    mistgenefilter.process(path, threshhold, testtype, testtypename, markerout)

def mistR(path, outpath, outfile, marker, testtypename):
    mistreptfilter.process(path, outpath, outfile, marker, testtypename)


def arguments():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subfunction')
    a_parser = subparsers.add_parser('Generate')
    a_parser.add_argument('-o', '--outpath', default='/home/cintiq/PycharmProjects/misty/mistout/', help='output folder for json files from mist')
    a_parser.add_argument('-t', '--testtype', required=True, help='path to test markers file')
    a_parser.add_argument('-T', '--testtypename', required=True, help='name of test, ex. CGF119')
    a_parser.add_argument('--threshhold', type=int, required=True, help='maximum amount of missed genes tolerated to allow a strain to pass')
    a_parser.add_argument('--fastaoutpath', default='./mistpass/', help='folder to place fasta file symlinks for those that pass the threshold')
    a_parser.add_argument('--distoutpath', default='./mistreport/', help='output folder for the report summary')
    a_parser.add_argument('--distoutfile', default='report', help='name of the report file')
    a_parser.add_argument('-a', '--alleles', default='/home/cintiq/Desktop/campylobacterjejuni/alleles/', help='folder that contains all allele files for mist requirements')
    a_parser.add_argument('-c', '--cores', default=multiprocessing.cpu_count(), help='number of cores to run')
    a_parser.add_argument('path', nargs='+', help='fasta files to run mist on')

    b_parser = subparsers.add_parser('Refine')
    b_parser.add_argument('--genomethreshhold', type=int, required=True, help='cutoff number for genomes that gene is not in')
    b_parser.add_argument('-t', '--testtype', required=True, help='path to and name of test markers file eg. MLST.markers')
    b_parser.add_argument('-T', '--testtypename', required=True, help='name of test to be performed eg. MLST')
    b_parser.add_argument('--markerout', default='/home/cintiq/PycharmProjects/misty/marker/', help='folder to place new marker in')
    b_parser.add_argument('--reptoutpath', default='./mistreport/', help='output folder for the report summary')
    b_parser.add_argument('--reptoutfile', default='refinedreport.json', help='name of the report file')
    b_parser.add_argument('-o', '--outpath', default='/home/cintiq/PycharmProjects/misty/mistout/', help='output folder for json files from mist')
    b_parser.add_argument('path', help='report file')
    return parser.parse_args()

def processGEN(path, outpath, testtype, testtypename, alleles, distoutpath, distoutfile, fastaoutpath, threshhold, cores):
    print('Running MIST...')
    mistM(path, outpath, testtypename, testtype, alleles, cores)
    print('Performing update...')
    jsonlist = glob.glob(outpath+'*./json')
    for json in jsonlist:
        mistU(alleles, json, testtypename)
    print('Making report...')
    letter = mistD(outpath, distoutpath, distoutfile, testtype, testtypename)
    distreploc=os.path.join(distoutpath, distoutfile+letter+'.json')
    print('Symlinking passing fasta files')
    mistF(distreploc, fastaoutpath, threshhold, outpath, testtypename)
    print('Generate has completed')


def processREF(path, genomethreshhold, testtype, testtypename, markerout, outpath, outfile):
    print('Filtering genes...')
    mistG(path, genomethreshhold, testtype, testtypename, markerout)
    print('Revising report based on new genes...')
    mistR(path, outpath, outfile, markerout, testtypename)
    print('Refine has completed')




def main():
    args = arguments()
    if args.subfunction=='Generate':
        processGEN(args.path, args.outpath, args.testtype, args.testtypename, args.alleles, args.distoutpath, args.distoutfile, args.fastaoutpath, args.threshhold, args.cores)
    if args.subfunction=='Refine':
        processREF(args.path, args.genomethreshhold, args.testtype, args.testtypename, args.markerout, args.reptoutpath, args.reptoutfile)




if __name__ == '__main__':
    main()