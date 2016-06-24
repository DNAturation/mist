# takes in the report file and outputs a csv file,
# giving numbers of remaining genes and genomes at certain numbers of threshold tolerance
import csv
import json
import os
import argparse

def pathfinder(outpath):
    '''makes output directory'''
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)

def reader(path):
    '''reads report.json file'''
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def gettotals(data):
    '''returns max amounts of genes and genomes missing, serves as an end point for analysis'''
    genomemax=max([len(data["GenomesMissingGenes"][x]) for x in data["GenomesMissingGenes"]])
    genemax=max([len(data["GenesMissingGenomes"][x]) for x in data["GenesMissingGenomes"]])

    print(genomemax, genemax)
    return genomemax, genemax



def genecull(genethresh, data):
    '''creates a dictionary of the current cutoff threshold and number of passing genes'''
    gened={}
    gened[genethresh]=0
    for genes in data['GenesMissingGenomes']: #iterate through gene names in the report file
        if len(data['GenesMissingGenomes'][genes]) <= genethresh:   #counts genes that pass (have fewer missing
            gened[genethresh] += 1                                  #genomes than the provided threshold)
    return gened


def genomecull(genomethresh, data):
    '''generates dictionary of the current cutoff threshold and number of passing genomes'''
    genomed={}
    genomed[genomethresh]=0
    for genomes in data['GenomesMissingGenes']: #iterate through genome names in the report file
        if len(data['GenomesMissingGenes'][genomes]) <= genomethresh:   #counts genomes that pastt(have fewer missing
            genomed[genomethresh] += 1                                  #genomes than the provided threshold)
    return genomed


def header(outpath, outfile):
    '''makes the csv file and writes the headers to it'''
    with open(os.path.join(outpath, outfile), 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        csvwriter.writerow(['GenomeThreshold', 'GenomesLeft', '', 'GeneThreshold', 'GenesLeft'])

def writer(outpath, outfile, genomeculld, geneculld):
    '''writes both genes and genomes, their cutoffs and number of passing genes/genomes'''
    with open(os.path.join(outpath, outfile), 'a', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        for (genomek, genomev), (genek, genev) in zip(sorted(genomeculld.items()), sorted(geneculld.items())): #iterates through both dictionaries
            csvwriter.writerow([genomek, genomev, '', genek, genev])

def writergenome(outpath, outfile, genomeculld):
    '''writes number of passing genomes and their cutoffs. Used when there are more genomes than genes(genes run out before genomes)'''
    with open(os.path.join(outpath, outfile), 'a', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        for (genomek, genomev) in sorted(genomeculld.items()): #iterates only through the genome dictionary
            csvwriter.writerow([genomek, genomev])

def writergene(outpath, outfile, geneculld):
    '''writes number of passing genomes and their cutoffs. Used when there are more genes than genomes(genomes run out before genes)'''
    with open(os.path.join(outpath, outfile), 'a', newline='') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        for (genek, genev) in sorted(geneculld.items()): #iterates only through the gene dictionary
            csvwriter.writerow(['', '', '', genek, genev])



def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genethreshhold', default=None, type=int, help='maximum number within a gene of missing genomes tolerated before gene is removed')
    parser.add_argument('--genomethreshhold', default=None, type=int, help='maximum number within a genome of missing genes tolerated before genome is removed')
    parser.add_argument('--genomemin', default=0, type=int, help='minimum number cutoff to start at for genome')
    parser.add_argument('--genemin', default=0, type=int, help='minimum number cutoff to start at for gene')
    parser.add_argument('--outfile', default='simulator.csv', help='name of file to be created')
    parser.add_argument('-o', '--outpath', default='./sim/', help='directory for the summary file')
    parser.add_argument('path', help='path to report file')
    return parser.parse_args()

def process(path, outpath, outfile, genemin, genethreshhold, genomemin, genomethreshhold):
    data = reader(path)
    pathfinder(outpath)

    genomegene = True   #flags for which dictionary to iterate over, necessary due to current structure
    genome=False        #iterating over both gene and genome dictionaries at the same time
    gene=False

    maxes = gettotals(data) #get maxes of genomes and genes missing to serve as a starting point for threshold cutoffs
    genomemax=maxes[0]
    genemax=maxes[1]
    if genethreshhold != None:
        if int(genethreshhold) > genemax:   #if input gene threshold is higher than the max missing genomes a gene has,
            genethreshhold=genemax          #set it to the max
    if genomethreshhold != None:
        if int(genomethreshhold) > genomemax:   #if input genome threshold is higher than the max missing genes a genome has,
            genomethreshhold=genomemax          #set it to the max

    currentcount=0 #counts number of iterations, and thus current threshold cutoff number
    leavingoff=0 #used to calculate where to begin again during switching to another writer due to genome or gene ending before the other

    header(outpath, outfile) #makes the csv file and writes the header to it

    while genomegene:
        for x, y in zip(range(int(genemin), int(genethreshhold)), range(int(genomemin), int(genomethreshhold))): #iterates through a range of the max of both thresholds
                                                                                                                #the second threshold may be unnecessary

            if x+1 == genomethreshhold and x+1 != genethreshhold: #reaches max missing of genomes but genes continue
                '''switch to genome writer'''
                geneculld=genecull(currentcount, data) #sets dict variable to a single length dictionary of key threshold
                genomeculld=genomecull(currentcount, data) #and value number of passing genes/genomes
                writer(outpath, outfile, genomeculld, geneculld) #writes to csv file from information in the dictionaries
                currentcount+=1 #increase counter, current threshold
                leavingoff+=currentcount #sets leavingoff to current threshold, for use in switching to another iterator
                gene=True #ets flag to iterate over gene dictionary
                genomegene=False #removes flag to iterate over both dictionaries
                break

            elif x+1 == genethreshhold and x+1 != genomethreshhold: #reaches max missing of genes but genomes continue
                '''switch to gene writer'''
                geneculld=genecull(currentcount, data) #same as above, but specifies genome flag instead of gene
                genomeculld=genomecull(currentcount, data)
                writer(outpath, outfile, genomeculld, geneculld)
                currentcount+=1
                leavingoff+=currentcount
                genome=True
                genomegene=False
                break

            elif x+1 == genethreshhold and x+1 == genomethreshhold:
                '''when both conclude at same time, end'''
                geneculld=genecull(currentcount, data)
                genomeculld=genomecull(currentcount, data)
                writer(outpath, outfile, genomeculld, geneculld)
                return

            else:
                '''an end has not been reached, continue the loop'''
                geneculld=genecull(currentcount, data)
                genomeculld=genomecull(currentcount, data)
                writer(outpath, outfile, genomeculld, geneculld)
                currentcount+=1

    while genome:
        for x in range(int(genomethreshhold) - leavingoff):#iterates through threshold, remaining after taking out what was done already
            genomeculld = genomecull(currentcount, data) #set genomeculld to length 1 dictionary
            writergenome(outpath, outfile, genomeculld) #write the dictionary to the csv file
            currentcount+=1 #increase the count
            if currentcount == genomethreshhold+1: #reaches end of genomethreshold, ends
                return

    while gene:
        for x in range(int(genethreshhold) - leavingoff):#iterates through threshold, remaining after taking out what was done already
            geneculld = genecull(currentcount, data) #set geneculld to length 1 dictionary
            writergene(outpath, outfile, geneculld) #write the dictionary to the csv file
            currentcount+=1 #increase the count
            if currentcount == genethreshhold+1: #reaches end of gene threshold, ends
                return

def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.genemin, args.genethreshhold, args.genomemin, args.genomethreshhold)


if __name__ == '__main__':
    main()