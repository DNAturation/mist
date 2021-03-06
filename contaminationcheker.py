# script written for checking how many genes are missing from a genome, after user determined cutoff, removes the genome
# and considers it as contamination
# Takes in report.json from MIST and returns a new json report file of MIST as well as a list of the
# contaminating genomes, also whether the strain is isolated from chicken or human
#
# note this is not for any pipeline, but was written for personal use for the term-end presentation;
# many things are hard-coded and don't make sense out of context

import argparse
import os
import json
import csv
import re
from Bio import SeqIO

def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)


def reader(path):
    with open(path, 'r') as f:
        data = json.load(f)

    return data


def checker(data, thresh):  # adds all genomes that have over threshold missing genes to a list; genomes that will be thrown out
    contaminatinggenomes=[]
    for strain in data['GenomesMissingGenes']:
        if len(data['GenomesMissingGenes'][strain]) >= int(thresh):
            contaminatinggenomes.append(strain)
    cont=['2015-SEQ-1679', '2015-SEQ-1680', '2015-SEQ-1682', '2013-SEQ-0010', '2013-SEQ-0011']
    return contaminatinggenomes+cont


def newjsonmaker(contaminatinggenomes, data):  # creates dictionary for new json MIST report
    newd = dict()
    newd['GenesMissingGenomes'] ={}
    newd['GenomesMissingGenes'] = {}
    for gene in data['GenesMissingGenomes']:
        newd['GenesMissingGenomes'][gene] = list()
        for genome in data['GenesMissingGenomes'][gene]:
            if genome not in contaminatinggenomes:
                newd['GenesMissingGenomes'][gene].append(genome)

    for genome in data['GenomesMissingGenes']:
        if genome not in contaminatinggenomes:
            newd['GenomesMissingGenes'][genome] = data['GenomesMissingGenes'][genome]
    return newd


def repwriter(contaminatinggenomes, newd, outpath, repname, conlist):
    with open(os.path.join(outpath, repname), 'w') as f:
        json.dump(newd, f)
    with open(os.path.join(outpath, conlist), 'w') as g:
        sep = ', '
        g.write(sep.join(contaminatinggenomes))


def remover(genomes, contaminatinggenomes):
    files = os.listdir(genomes)
    for file in files:
        if os.path.splitext(file)[0][:-6] in contaminatinggenomes:
            # print('removing file', file)
            os.remove(os.path.join(genomes, file))


def accessorybinaryremover(contaminatinggenomes, fapath, newfa):  # removes genomes from the roary files
    with open(fapath, 'r') as f:
        with open(newfa, 'w') as g:
            for record in SeqIO.parse(f, 'fasta'):
                if str(record.id) not in contaminatinggenomes:
                    SeqIO.write(record, g, 'fasta')
                else:
                    print('skipping', record.id)


def chickenorhuman(straindatapath):  # determines whether the strain came from humans or chickens based off the name
    chickenornot=dict()
    notchickenmatch = '^\d\d-'
    with open(straindatapath, 'r') as f:
        filereader = csv.reader(f, dialect='excel')
        for row in filereader:
            if re.match(notchickenmatch, row[1]):
                chickenornot[row[0]] = 1
            else:
                chickenornot[row[0]] = 0
    return chickenornot


def isolatedatamaker(data, contaminationlist, chickenornot):  # from all the genomes, where the stuff is from
    l=list()
    for genome in data['GenomesMissingGenes']:
        if genome not in contaminationlist:
            l.append([genome, chickenornot[genome]])
    return l


def csvwriter(l, isolatepath):
    header = ['Sequence identifier', 'is chicken']
    with open(isolatepath, 'w')as f:
        filewriter = csv.writer(f, delimiter='\t')
        filewriter.writerow(header)
        for listpair in l:
            filewriter.writerow(listpair)


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='MIST report json file to check')
    parser.add_argument('-t', '--thresh', default=7, help='how many missing genes before genome is considered contamination')
    parser.add_argument('--outpath', default='./newreport', help='path to dump all output files')
    parser.add_argument('--repname', default='newmistreport.json', help='name of the new json report file')
    parser.add_argument('--conlist', default='conlist.txt', help='name of file to hold contaminating genome names')
    parser.add_argument('--genomes', default='/home/phac/kye/autocreate/edassemblies/')
    parser.add_argument('--fapath', default='/home/phac/kye/assemblies_for_ed/results/pangenome_1468951604/accessory_binary_genes.fa')
    parser.add_argument('--newfa', default='/home/phac/kye/presentationstuff/edassem/accessory_binary_genes2.fa')
    parser.add_argument('--isolatekey', default='/home/phac/kye/presentationstuff/edassem/isolatekey.csv')
    parser.add_argument('--isolatepath', default='/home/phac/kye/presentationstuff/edassem/isolatedata.csv')
    return parser.parse_args()

def main():
    args = arguments()
    pathfinder(args.outpath),
    data = reader(args.path)
    contaminatinggenomes = checker(data, args.thresh)
    newd = newjsonmaker(contaminatinggenomes, data)
    repwriter(contaminatinggenomes, newd, args.outpath, args.repname, args.conlist)
    remover(args.genomes, contaminatinggenomes)
    accessorybinaryremover(contaminatinggenomes, args.fapath, args.newfa)
    chickenornot = chickenorhuman(args.isolatekey)
    chickend = isolatedatamaker(data, contaminatinggenomes, chickenornot)
    csvwriter(chickend, args.isolatepath)

if __name__ == '__main__':
    main()