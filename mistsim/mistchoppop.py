#creates chopped.csv files, which contains which genes remain after removing the x worst genes
import csv
import json
import os
import argparse
import multiprocessing

def pathfinder(outpath, chopout, startchop, rankedlist):
    '''
    creates output directory for symlinked genes after the chopping
    '''
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)
    if startchop:
        if not os.access(os.path.join(outpath, chopout), os.F_OK):
            os.mkdir(os.path.join(outpath, chopout))
        for start in startchop:
            if int(start) < len(rankedlist):
                if not os.access(os.path.join(outpath, chopout, start)+'/', os.F_OK):
                    os.mkdir(os.path.join(outpath, chopout, start)+'/')
            else:
                break
    else: #if no start point is given, defaults to creating just a 0, containing all genes
        if not os.access(os.path.join(outpath, chopout), os.F_OK):
            os.mkdir(os.path.join(outpath, chopout))
        if not os.access(os.path.join(outpath, chopout, '0'), os.F_OK):
            os.mkdir(os.path.join(outpath, chopout, '0'))

def jreader(path):
    '''
    opens report.json file, used for the basis of looking up missing genes from genomes and how bad a gene is
    '''
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def genemax(data):
    '''
    finds the maximum number of genomes missing from a gene, serves as start point for the threshold cut down
    '''
    genemax=max([len(data["GenesMissingGenomes"][x]) for x in data["GenesMissingGenomes"]])
    return genemax

def ranker(data, gmax, startpop, endpop):
    '''
    organizes the genes present into a descending list based on the number of genes they are missing
    (most missing are the first ones)
    '''
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


def chopper(rankedlist, startchop):
    '''
    takes in the list of all genes (ranked from worst to best), and returns which genes are left given
    how many genes to remove
    '''
    if startchop:
        chopped = rankedlist[int(startchop):]
        return chopped
    else:
        chopped = rankedlist[startchop:]
        return chopped

def symmer(chopped, outpath, chopout, alleles, cores, start):
    '''
    does the symlinking of files in parallel
    '''
    if start:
        pool = multiprocessing.Pool(int(cores))
        for gene in chopped:
            if gene + '.fasta' not in os.listdir(os.path.join(outpath, chopout, start)+'/'):
                pool.apply_async(os.symlink, args=(os.path.join(alleles, gene)+'.fasta',
                                            os.path.join(outpath, chopout, start)+'/' + gene +'.fasta'))
        pool.close()
        pool.join()
    else: #if no start argument is given, uses 0 instead of None to prevent errors
        pool = multiprocessing.Pool(int(cores))
        for gene in chopped:
            if gene + '.fasta' not in os.listdir(os.path.join(outpath, chopout, '0')):
                pool.apply_async(os.symlink, args=(os.path.join(alleles, gene)+'.fasta',
                                            os.path.join(outpath, chopout, '0', gene)+'.fasta'))
        pool.close()
        pool.join()


def writer(outpath, outfile, chopped, start):
    '''
    writes out the reamining genes to a csv file
    '''
    x=str(start)+outfile
    with open(os.path.join(outpath, x), 'w') as f:
        csvwriter = csv.writer(f, delimiter=' ')
        for gene in chopped:
            gene = gene.split()
            csvwriter.writerow(gene)

def arguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('--outpath', default='./sim/')
    parser.add_argument('--outfile', default='chopped.csv')
    parser.add_argument('--startpop', default=None)
    parser.add_argument('--endpop', default=None)
    parser.add_argument('--startchop', default=None, nargs='+')
    parser.add_argument('--endchop', default=None)
    parser.add_argument('--chopout', default='chopsyms/')
    parser.add_argument('--alleles', default='~/kye/autocreate/prissy/alleles/')
    parser.add_argument('--cores', default=int(multiprocessing.cpu_count()))
    parser.add_argument('path', help='report json file')
    return parser.parse_args()


def process(path, outpath, outfile, startpop, endpop, startchop, chopout, alleles, cores):
    data = jreader(path)
    gmax = genemax(data)
    rankedlist= ranker(data, gmax, startpop, endpop)
    pathfinder(outpath, chopout, startchop, rankedlist)
    if startchop:
        for start in startchop: #startchop is an nargs = + argument, so it's expecting a list with multiple start points
            if int(start) < len(rankedlist):#incase the inputted starting points exceed the number of genes,
                                            # this keeps it from looping more than necessary
                if not os.path.isfile(os.path.join(outpath, (str(start)+outfile))): #prevents writing a new gene
                                                                                    # file if one already exists
                    chopped = chopper(rankedlist, start)
                    symmer(chopped, outpath, chopout, alleles, cores, start)
                    writer(outpath, outfile, chopped, start)
                else:
                    print('Skipping ', str(start), 'due to file existing')
            else:
                break

    else:
        chopped = chopper(rankedlist, startchop)
        symmer(chopped, outpath, chopout, alleles, cores, startchop)
        writer(outpath, outfile, chopped, startchop)

def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.startpop, args.endpop, args.startchop, args.chopout, args.alleles, args.cores)

if __name__ == '__main__':
    main()