#pop worst genes off one at a time and look to see if the genomes looks better (perfect)
import csv
import json
import os
import argparse

def pathfinder(outpath):
    '''creates output directory'''
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)


def jreader(path):
    '''opens report.json file, used for the basis of looking up missing genes from genomes and how bad a gene is'''
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def genemax(data):
    '''finds the maximum number of genomes missing from a gene, serves as start point for the threshold cut down'''
    genemax=max([len(data["GenesMissingGenomes"][x]) for x in data["GenesMissingGenomes"]])
    return genemax

def setmaker(data):
    '''makes a dictionary of key = genomes, value = set of missing genes, from which the values will be popped
    in order to obtain pure genomes'''
    dgenome = {}
    for genome in data["GenomesMissingGenes"]: #access the names of genomes in the report.json file
        dgenome[genome]=set(data["GenomesMissingGenes"][genome])    #create new dictionary of key genome name,
                                                                    #value list of missing genes for use in popping
                                                                    #out genes in that list to find perfect genomes
    return dgenome

def ranker(data, gmax, startpop, endpop):
    '''organizes the genes present into a descending list based on the number of genes they are missing (most missing
    are the first ones)'''
    rankedlist = []
    for thresh in range(gmax+1)[startpop:endpop:-1]:#iterates through a number(thresh = threshold) from highest (default gmax)
                                                    # to lowest (Default 0) in order to determine which genes to take out

        for genes in data["GenesMissingGenomes"]: #iterates through the names of the genes from the report.json file

            if len(data["GenesMissingGenomes"][genes]) == thresh:   #matches the number of genomes missing for the gene to the
                rankedlist.append(genes)                            #thresh, which starts at the highest number. The genes
                                                                    #missing from the most genomes get chosen and appended to
                                                                    #a list of genes first, resulting in a list that starts with
                                                                    #the worst genes and ends in the best genes

    return rankedlist

def writelister(rankedlist, dgenome):
    '''makes a list, with the number of genes removed and the number of genomes without missing genes, to be written
    to csv file'''
    writelist = []
    for index, gene in enumerate(rankedlist): #iterates through the ordered list of the worst genes
        puregenomes=0   #resets 'puregenomes' count every iteration of a new gene
        for genome in dgenome:  #looks through the made dictionary of key genomes value list of genes and removes
                                #genes in the genome's list that match the current iteration
            try:
                dgenome[genome].remove(gene)
            except:
                pass

            if dgenome[genome] == set():#checks for genomes that have no missing genes, if found, adds one for each
                puregenomes+=1          #present genome to the puregenomes counter

        writelist.append([index, puregenomes])  #makes a list of two items, number of genes removed and number of genomes
                                                #that are pure for future csv file writing

    return writelist


def writer(outpath, outfile, writelist):
    '''writes to csv file the number of genes removed and the number of genomes that contain all remaining genes'''
    with open(os.path.join(outpath, outfile), 'a') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        csvwriter.writerow(writelist)

def header(outpath, outfile):
    with open(os.path.join(outpath, outfile), 'w') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        csvwriter.writerow(['NumberOfGenesRemoved', 'NumberOfPerfectGenomes'])


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outpath', default='./sim/')
    parser.add_argument('--outfile', default='pop.csv')
    parser.add_argument('--startpop', default=None)
    parser.add_argument('--endpop', default=None)
    parser.add_argument('path')
    return parser.parse_args()


def process(path, outpath, outfile, startpop, endpop):
    data = jreader(path) #reads report.json file and saves into data
    pathfinder(outpath) #make output directory

    gmax = genemax(data)#obtains the maximum number of missing genomes from genes
    if startpop != None:
        if startpop > gmax: #if entered starting threshold is greater than the maximum number of missing genomes from a
            startpop = gmax # gene, sets it to the maximum number of missing genomes

    dgenome = setmaker(data) #stores temporary dictionary to pop things out of in the dgenome variable
    rankedlist=ranker(data, gmax, startpop, endpop) #ordered list of worst genes
    writelist = writelister(rankedlist, dgenome) #creates a 2D list for what to write to csv file
    header(outpath, outfile) #makes the initial popper csv file and writes the headers
    for lines in writelist: #iterate through each line to write in the list of lines to write and writes them
        writer(outpath, outfile, lines)


def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.startpop, args.endpop)


if __name__ == '__main__':
    main()