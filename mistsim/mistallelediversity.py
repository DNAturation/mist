#made to see which alleles of each gene are represented the most; makes a json file with key = gene, value = dictionary
#, and from that dictionary, key = allele, value = number of times it shows up

import argparse
import csv
import os
import glob
import json
from Bio import SeqIO


def pathfinder(output):
    if not os.access(output, os.F_OK):
        os.mkdir(output)

def genegetter(outpath, file):
    '''
    gets the names of all the genes in order of worst to best
    '''
    with open(outpath + file + 'chopped.csv', 'r') as f: #this is a csv file with a single column of names
        namefile = csv.reader(f, delimiter=' ')
        genenames=[]
        for line in namefile:
            genenames.append(line)
    return genenames

def allelegetter(outpath, chopout, file, genenames, start, end):
    '''
    creates a dictionary of every allele in each gene
    '''
    dalleles={}
    for name in genenames[int(start):int(end)+1]: #picks out which genes to actually use
        name = name[0] #the name gets returned as a list, this is to turn it into a string for further use
        dalleles[name]={}
        for allele in SeqIO.parse(os.path.join(outpath, chopout, file, name) +'.fasta', 'fasta'): #access each allele in the gene's fasta file
            dalleles[name][allele.id]=allele.seq
    return dalleles

def genomegetter(pristine):
    '''
    gets a list of genomes to test on
    '''
    genomes = os.listdir(pristine)
    return genomes

def reportreader(reportdir, genomes):
    '''
    reads json file for every genome
    '''
    for genome in genomes:
        report = reportdir+genome
        with open(report, 'r') as f:
            jsonfile=json.load(f)
        yield jsonfile


def counter(jsonfile, genenames, dalleles, start, end):
    '''
    counts the number of occurrances of an allele within genomes
    '''
    ddiv={}
    for dict in jsonfile: #iterates through every report json file for each strain
        for name in genenames[int(start):int(end)+1]: #picks out the genes to actually use
            name = name[0] #name gets returned as a length 1 list, this converts it to a string
            if name not in ddiv.keys():
                ddiv[name]={}
            for id in dalleles[name]: #access the id number for each allele for each gene
                # print(dict["Results"][0]["TestResults"]["wgmlst"][name]["BlastResults"]["QueryAln"])
                if dict["Results"][0]["TestResults"]["wgmlst"][name]["BlastResults"]["QueryAln"] == dalleles[name][str(id)]: #sequence for the associated allele number
                    try:
                        ddiv[name][str(id)]+=1
                    except:
                        ddiv[name][str(id)]=1

    return ddiv


def writer(ddiv, output, outfile, start, end):
    with open(os.path.join(output, str(start)+str(end)+outfile), 'w+') as f:
        json.dump(ddiv, f, indent=4, sort_keys=True)

def arguments():
    parser= argparse.ArgumentParser()
    parser.add_argument('--file', default='0', help='lowest chopped.csv file, containing as many genes as needed')
    parser.add_argument('--chopout', default='chopsyms/', help='access the directory of symlinked gene files')
    parser.add_argument('--reportdir', default='/home/phac/kye/misty/mistout/', help='directory of all genome report files generated by MIST')
    parser.add_argument('--pristinedir', default='/home/phac/genomes/campylobacter/pristine/', help='pristine genomes directory, which this is made to test')
    parser.add_argument('--output', default='diversity/', help='where to place the files')
    parser.add_argument('--outfile', default='diversity.json', help='name of the output json file')
    parser.add_argument('--start', default=0, help='where to start, in terms of worst gene to compare')
    parser.add_argument('--end', default=10, help='where to end the comparison of gene alleles')
    return parser.parse_args()

def process(output, outfile, outpath, file, chopout, reportdir, pristine, start, end):
    pathfinder(output)
    print('grabbing gene names')
    genenames = genegetter(outpath, file)
    print('retrieving alleles')
    dalleles=allelegetter(outpath, chopout, file, genenames, start, end)
    print('retrieving genome names')
    genomes = genomegetter(pristine)
    jsonfile = reportreader(reportdir, genomes)
    print('counting...')
    ddiv = counter(jsonfile, genenames, dalleles, start, end)
    print('writing to file')
    writer(ddiv, output, outfile, start, end)

def main():
    args = arguments()
    process(args.output, args.outfile, args.outpath, args.file, args.chopout, args.genomes,
            args.pristine, args.start, args.end)

if __name__ == '__main__':
    main()
