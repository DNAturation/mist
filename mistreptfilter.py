# removes genes and species from the results file that are no longer in the new markers cutoff file
import json
import os
import argparse
import re
import csv

def pathfinder(outpath):
    '''makes output dir'''
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)

def reader(path):
    '''opens the report file'''
    with open(path, 'r') as f:
        data = json.load(f)
    return data


def glist(marker, testtypename):
    '''opens the new markers file and retrieves all the gene names'''
    with open(os.path.join(marker, testtypename)+'.markers', 'r') as f:
        try:  # check to see if the marker is in json format
            genelist = []
            for key in json.load(f).keys():  # iterates through all gene names in the markers file
                genelist.append(key)  # collects all gene names in genelist
        except ValueError:
            f.seek(0)
            genelist=[]
            data = f.readlines()[1:]  # skips the first header line and reads the data
            data = [line.split()[0] for line in data]  # gets the names only of genes from new markers file
            if data != []:  # skips empty names
                for line in data:  # takes the name out from the list
                    genelist.append(line)  # adds names to genelist
    return genelist

def gatherer(data, genelist):
    '''creates the new set of genes and genomes'''
    d = {}
    GenesMissingGenomes={}
    GenomesMissingGenes={}
    for gene in data["GenesMissingGenomes"]:
        if gene in genelist: #checks if the name of the gene is in the new markers file, if so, adds it to the dictionary for the gene portion
            GenesMissingGenomes[gene] = data["GenesMissingGenomes"][gene]
    for genome in data["GenomesMissingGenes"]:
        GenomesMissingGenes[genome] = []
        for mgenes in data["GenomesMissingGenes"][genome]:
            if mgenes in genelist:
                GenomesMissingGenes[genome].append(mgenes) #adds every gene present in the markers file to a particular genome
    d['GenesMissingGenomes'] = GenesMissingGenomes#addes the two dictionaries into one larger dictionary for writing
    d['GenomesMissingGenes'] = GenomesMissingGenes
    return d



def writer(outpath, outfile, d):
    with open(os.path.join(outpath, outfile), 'a') as f:
        json.dump(d, f, indent=4, sort_keys=True)

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outpath', default='/home/cintiq/PycharmProjects/misty/mistreport/', help='output directory for the report')
    parser.add_argument('-n', '--outfile', default='refinedreport.json', help='output file name')
    parser.add_argument('--marker', default='/home/cintiq/PycharmProjects/misty/marker/', help='directory containing new .markers file of interest')
    parser.add_argument('--testtypename', required=True, help='the name of the test run')
    parser.add_argument('path')
    return parser.parse_args()

def process(path, outpath, outfile, marker, testtypename):
    pathfinder(outpath)
    data = reader(path)
    genelist = glist(marker, testtypename)
    d = gatherer(data, genelist)
    writer(outpath, outfile, d)

def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.marker, args.testtypename)

if __name__ == '__main__':
    main()