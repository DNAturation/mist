#reads report generated by mistdistrept and looks at the genes that don't aren't present in genomes. if genes aren't in genomes over a threshold(defined by user), the gene is not copied to a new temporary marker file

import argparse
import os
import json

def pathfinder(markers):
    if not os.access(markers, os.F_OK):
        os.mkdir(markers)

def reader(path):
    '''opens the report file made by mistdistrept and stores it into data variable'''
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def countergenomes(data):
    genelist = data['GenesMissingGenomes']
    for gene in genelist:
        missinggno = len(genelist[gene])
        yield missinggno, gene

def filter(missinggno, gene, threshhold, blacklist):
    '''determines which genes fail the cutoff, and appends them to a list that is created outside
    this function within the process function'''
    if missinggno > threshhold:
        blacklist.append(gene)

def cull(blacklist, testtype, testtypename, markers):
    '''writes the new marker data out if they are not in the list of removed genes/genomes'''
    with open(testtype, 'r') as ref:
        with open(os.path.join(markers, testtypename)+'.markers', 'a+') as temp:
            lines = ref.readlines()
            for line in lines:
                if blacklist == []:
                    temp.write(line)
                if not any(gene in line for gene in blacklist):
                    temp.write(line)




def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-thresh', '--threshhold', type=int, required=True, help='cutoff number for genomes that gene is not in')
    parser.add_argument('-T', '--testtypename', required=True, help='type of test performed eg. MLST')
    parser.add_argument('-t', '--testtype', required=True, help='path to and name of test markers file')
    parser.add_argument('-markers', '--markersout', default='/home/cintiq/PycharmProjects/misty/markers/', help='markers folder to place new .marker in')
    parser.add_argument('path', help='report file')
    return parser.parse_args()

def process(path, threshhold, testtype, testtypename, markers):
    blacklist = []
    pathfinder(markers)
    data = reader(path)
    x = countergenomes(data)
    for missinggno, gene in x:
        filter(missinggno, gene, threshhold, blacklist)
    cull(blacklist, testtype, testtypename, markers)

def main():
    args = arguments()
    process(args.path, args.threshhold, args.testtype, args.testtypename, args.markerspout)

if __name__ == '__main__':
    main()