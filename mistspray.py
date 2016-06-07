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

def mistM(path, outpath, testtype, mistfolder):
    '''runs mist'''
    mistmain.process(path, outpath, testtype, mistfolder)

def mistF(path, fastaoutpath, threshhold, outpath, testtype):
    '''runs a userdefined threshold cutoff and symlinks fasta files that pass'''
    mistfastasym.process(path, fastaoutpath, threshhold, outpath, testtype)

def mistD(path, outpath, outfile, testtype, mistfolder):
    '''generates a report'''
    mistdistrept.process(path, outpath, outfile, testtype, mistfolder)

def mistU(alleles, mistfolder, testtype):
    '''runs updater on mist main output, dillon's script'''
    update_definitions.process(alleles, mistfolder, testtype)

def mistG(path, threshhold, mistfolder, testtype, markerout):
    '''makes a temporary .markers file based on genes that fall under a threshold of not being present in genomes'''
    mistgenefilter.process(path, threshhold, mistfolder, testtype, markerout)

def mistR(path, outpath, outfile, marker, testtype):
    mistreptfilter.process(path, outpath, outfile, marker, testtype)


def arguments():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subfunction')
    a_parser = subparsers.add_parser('Generate')
    a_parser.add_argument('-o', '--outpath', default='/home/cintiq/PycharmProjects/misty/mistout/', help='output folder for json files from mist')
    a_parser.add_argument('-m', '--mistfolder', default='/home/cintiq/Desktop/campylobacterjejuni/', help='folder that contains all mist marker file requirements')
    a_parser.add_argument('-t', '--testtype', required=True, help='type of test, ex. CGF119')
    a_parser.add_argument('--threshhold', type=int, required=True, help='maximum amount of missed genes tolerated to allow a strain to pass')
    a_parser.add_argument('--fastaoutpath', default='./mistpass/', help='folder to place fasta file symlinks for those that pass the threshold')
    a_parser.add_argument('--distoutpath', default='./mistreport/', help='output folder for the report summary')
    a_parser.add_argument('--distoutfile', default='report.json', help='name of the report file')
    a_parser.add_argument('-a', '--alleles', default='/home/cintiq/Desktop/campylobacterjejuni/alleles/', help='folder that contains all allele files for mist requirements')
    a_parser.add_argument('path', nargs='+', help='fasta files to run mist on')

    b_parser = subparsers.add_parser('Refine')
    b_parser.add_argument('--genomethreshhold', type=int, required=True, help='cutoff number for genomes that gene is not in')
    b_parser.add_argument('-t', '--testtype', required=True, help='type of test performed eg. MLST')
    b_parser.add_argument('--markerout', default='/home/cintiq/PycharmProjects/misty/marker/', help='folder to place new marker in')
    b_parser.add_argument('--reptoutpath', default='./mistreport/', help='output folder for the report summary')
    b_parser.add_argument('--reptoutfile', default='refinedreport.json', help='name of the report file')
    b_parser.add_argument('-o', '--outpath', default='/home/cintiq/PycharmProjects/misty/mistout/', help='output folder for json files from mist')
    b_parser.add_argument('--mistfolder', default='/home/cintiq/Desktop/campylobacterjejuni/', help='folder that contains all mist marker file requirements')
    b_parser.add_argument('path', help='report file')
    return parser.parse_args()

def processGEN(path, outpath, testtype, mistfolder, alleles, distoutpath, distoutfile, fastaoutpath, threshhold):
    mistM(path, outpath, testtype, mistfolder)
    mistU(alleles, outpath, testtype)
    mistD(outpath, distoutpath, distoutfile, testtype, mistfolder)
    distreploc=os.path.join(distoutpath,distoutfile)
    mistF(distreploc, fastaoutpath, threshhold, outpath, testtype)


def processREF(path, genomethreshhold, mistfolder, testtype, markerout, outpath, outfile):
    mistG(path, genomethreshhold, mistfolder, testtype, markerout)
    mistR(path, outpath, outfile, markerout, testtype)




def main():
    args = arguments()
    if args.subfunction=='Generate':
        processGEN(args.path, args.outpath, args.testtype, args.mistfolder, args.alleles, args.distoutpath, args.distoutfile, args.fastaoutpath, args.threshhold)
    if args.subfunction=='Refine':
        processREF(args.path, args.genomethreshhold, args.mistfolder, args.testtype, args.markerout, args.reptoutpath, args.reptoutfile)
#     # generate
# # if args.Generate == 'Generate':
#     mistM(args.path, args.outpath, args.testtype, args.mistfolder)
#     mistU(args.alleles, args.outpath, args.testtype)
#     mistD(args.outpath, args.distoutpath, args.distoutfile, args.testtype, args.mistfolder)
#     distreploc=os.path.join(args.distoutpath,args.distoutfile)
#     mistF(distreploc, args.fastaoutpath, args.threshhold, args.outpath, args.testtype)
#
#     # refine
# # if args.Refine == 'Refine':
#     mistG(args.path, args.genomethreshhold, args.mistfolder, args.testtype, args.markerout)
#     mistD(args.outpath, args.genomedistoutpath, args.genomedistoutfile, args.testtype, args.markerout)



if __name__ == '__main__':
    main()