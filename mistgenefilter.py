#reads report generated by mistdistrept and looks at the genes that don't aren't present in genomes. if genes aren't in genomes over a threshold(defined by user), the gene is not copied to a new temporary marker file

import argparse
import os
import json
import csv
import re

def pathfinder(markers):
    if not os.access(markers, os.F_OK):
        os.mkdir(markers)

def reader(path):
    '''opens the report file made by mistdistrept and stores it into data variable'''
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def testnamegetter(testtype):
    with open(testtype, 'r') as f:
        try: #try accessing markers file as .json file first
            data = json.load(f)
            for genome, keys in data.items():
                for key in keys:
                    if re.match('T(est)?\.?[-\._ ]?Name.*', key, flags=re.IGNORECASE):
                        return keys[key]

        except KeyError: #if access as .json file fails, try to access as csv file
            f.seek(0)
            reader=csv.reader(f, delimiter='\t')
            next(reader, None)
            for x in reader:
                testname=x[1]
                return testname

def countergenomes(data):
    genelist = data['GenesMissingGenomes'] #gets a list of every gene in the report file and their missing genomes
    for gene in genelist:
        missinggno = len(genelist[gene])    #gets the number of genomes a gene is missing from,
        yield missinggno, gene              #returns it along with the gene name

def filter(missinggno, gene, threshhold, blacklist):
    '''
    determines which genes fail the cutoff, and appends them to a list that is created outside
    this function within the process function
    '''
    if missinggno > threshhold:
        blacklist.append(gene)

def cull(blacklist, testtype, testtypename, markers):
    '''writes the new marker data out if they are not in the list of removed genes/genomes'''
    with open(testtype, 'r') as ref:  # open the original markers file
        with open(os.path.join(markers, testtypename)+'.markers', 'w') as temp:  # opens the new markers file for writing
            try:  # try loading the markers file as a json
                file = json.load(ref)
                tempd = {}
                for gene, geneinfo in file.items():  # loops through the gene name and other marker information
                    if gene not in blacklist:  # check to see if the gene passed, if so then adds it to the dictionary
                        tempd[gene] = geneinfo
                json.dump(tempd, temp)

            except ValueError: # use original csv style if .json form fails
                ref.seek(0)
                lines = ref.readlines() #prepares to write every line in original markers file to new markers file
                for line in lines:
                    if blacklist == []: #if the blacklist contains no genes, write every line
                        temp.write(line)
                    elif not any(gene in line for gene in blacklist): #compares the name of the genes present in blacklist
                        temp.write(line)                              #to the line to write, if not present in blacklist,
                                                                      #write the line


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-thresh', '--threshhold', type=int, required=True, help='cutoff number for genomes that gene is not in')
    parser.add_argument('-t', '--testtype', required=True, help='path to and name of test markers file')
    parser.add_argument('-markers', '--markersout', default='/home/cintiq/PycharmProjects/misty/markers/', help='markers folder to place new .marker in')
    parser.add_argument('path', help='report file')
    return parser.parse_args()

def process(path, threshhold, testtype, markers):
    blacklist = []
    pathfinder(markers)
    data = reader(path)
    testtypename = testnamegetter(testtype)
    x = countergenomes(data)
    for missinggno, gene in x:
        filter(missinggno, gene, threshhold, blacklist)
    cull(blacklist, testtype, testtypename, markers)

def main():
    args = arguments()
    process(args.path, args.threshhold, args.testtype, args.markerspout)

if __name__ == '__main__':
    main()