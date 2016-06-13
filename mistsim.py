# takes in the report file and outputs a csv file,
# giving numbers of remaining genes and genomes at certain numbers of threshold tolerance
import csv
import json
import os
import argparse

def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)

def reader(path):
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def genecull(genethresh, data):
    gened={}
    gened[genethresh]=0
    for genes in data['GenesMissingGenomes']:
        if len(data['GenesMissingGenomes'][genes]) <= genethresh:
            gened[genethresh] += 1
    return gened


def genomecull(genomethresh, data):
    '''generates dictionary of genomes '''
    genomed={}
    genomed[genomethresh]=0
    for genomes in data['GenomesMissingGenes']:
        if len(data['GenomesMissingGenes'][genomes]) <= genomethresh:
            genomed[genomethresh] += 1
    return genomed


def header(outpath, outfile):
    with open(os.path.join(outpath, outfile), 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        csvwriter.writerow(['GenomeThreshold', 'GenomesLeft', '', 'GeneThreshold', 'GenesLeft'])

def writer(outpath, outfile, genomeculld, geneculld):
    with open(os.path.join(outpath, outfile), 'a', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        for (genomek, genomev), (genek, genev) in zip(sorted(genomeculld.items()), sorted(geneculld.items())):
            csvwriter.writerow([genomek, genomev, '', genek, genev])

def writergenome(outpath, outfile, genomeculld):
    with open(os.path.join(outpath, outfile), 'a', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        for (genomek, genomev) in sorted(genomeculld.items()):
            csvwriter.writerow([genomek, genomev])

def writergene(outpath, outfile, geneculld):
    with open(os.path.join(outpath, outfile), 'a', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        for (genek, genev) in sorted(geneculld.items()):
            csvwriter.writerow(['', '', '', genek, genev])



def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genethreshhold', default=2000, help='maximum number within a gene of missing genomes tolerated before gene is removed')
    parser.add_argument('--genomethreshhold', default=40, help='maximum number within a genome of missing genes tolerated before genome is removed')
    parser.add_argument('--outfile', default='simulator.csv', help='name of file to be created')
    parser.add_argument('-o', '--outpath', default='./sim/', help='directory for the summary file')
    parser.add_argument('path', help='path to report file')
    return parser.parse_args()

def process(path, outpath, outfile, genethreshhold, genomethreshhold):
    data = reader(path)
    pathfinder(outpath)
    genomegene = True
    genome=False
    gene=False
    currentcount=1
    leavingoff=0
    header(outpath, outfile)
    while genomegene:
        for x, y in zip(range(genethreshhold), range(genomethreshhold)):
            if x+1 == genomethreshhold and x+1 != genethreshhold:
                geneculld=genecull(currentcount, data)
                genomeculld=genomecull(currentcount, data)
                writer(outpath, outfile, genomeculld, geneculld)
                currentcount+=1
                leavingoff+=currentcount
                gene=True
                genomegene=False
                break
            elif x+1 == genethreshhold and x+1 != genomethreshhold:
                geneculld=genecull(currentcount, data)
                genomeculld=genomecull(currentcount, data)
                writer(outpath, outfile, genomeculld, geneculld)
                currentcount+=1
                leavingoff+=currentcount
                genome=True
                genomegene=False
                break
            elif x+1 == genethreshhold and x+1 == genomethreshhold:
                geneculld=genecull(currentcount, data)
                genomeculld=genomecull(currentcount, data)
                writer(outpath, outfile, genomeculld, geneculld)
                return
            else:
                geneculld=genecull(currentcount, data)
                genomeculld=genomecull(currentcount, data)
                writer(outpath, outfile, genomeculld, geneculld)
                currentcount+=1

    while genome:
        for x in range(int(genomethreshhold) - leavingoff):
            genomeculld = genomecull(currentcount, data)
            writergenome(outpath, outfile, genomeculld)
            currentcount+=1
            if currentcount == genomethreshhold+1:
                return
    while gene:
        for x in range(int(genethreshhold) - leavingoff):
            geneculld = genecull(currentcount, data)
            writergene(outpath, outfile, geneculld)
            currentcount+=1
            if currentcount == genethreshhold+1:
                return

def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.genethreshhold, args.genomethreshhold)


if __name__ == '__main__':
    main()