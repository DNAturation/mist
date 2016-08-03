# takes in a directory of MIST output files (jsons) and creates a distance calculation of the alleles.
# If two strains have similar alleles, their numbers (intersection point on CSV file) are lower. Counts similarity
# by checking which allele MIST matched each read with; different alleles increase the score
# Note that this script assigns 'NA' values when a match for a gene could not be found, and these values
# contribute to the score (matching with other NA's). Also groups strains depending on their allele matches
import json
import argparse
import os
import glob
import subprocess
import collections
import tempfile
import csv


def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)


def filegetter(path):
    files = glob.glob(path + "*.json")
    return files


def reader(file, genedic, genedic2):
    # generates two dictionarys of the allele number that a strain's gene was matched with.
    # Note: contig truncations and no matches are counted as an NA allele
    with open(file, 'r') as f:
        d = {}
        data = json.load(f)
        strainname = os.path.splitext(os.path.basename(file))[0]
        results = data['Results'][0]['TestResults'] # access base important level in mist json output
        for test in results: # testname, used for accessing further values (eg. wgmlst)
            for gene in results[test]:
                if results[test][gene]["BlastResults"] is None or results[test][gene]["IsContigTruncation"]:
                    # treats contig truncations and no matches as an 'NA' allele match. Note that this will count
                    # towards the score
                    d[gene] = "NA"
                else:
                    # grabs which allele MIST matched to
                    d[gene] = results[test][gene]["AlleleMatch"]
                genedic[strainname] = d # updates the dictionaries (that are outside this function)
                genedic2[strainname] = d


def distget(gene1, gene2, dist):
    #increases distance score on difference in allele match
    if gene1 != gene2:
        dist += 1
    return dist


def compare(genedic, genedic2):
    # dmat is a dictionary that will be turned into a matrix.
    # It holds key of strain1, value of a dictionary(key = strain2, value = score)
    dmat = collections.defaultdict(dict)
    for strain1 in genedic:
        genedic2.pop(strain1)
        for strain2 in genedic2:
            dist = 0
            if strain1 != strain2:
                for gene1 in genedic[strain1]:
                    for gene2 in genedic[strain2]:
                        dist = distget(genedic[strain1][gene1], genedic2[strain2][gene2], dist)
                        dmat[strain1][strain2] = dist
    return dmat


def csvwriter(dmat, outpath, outfile):
    # outsources the prettyifying of the dictionary into a table to an R script, alleledistance.R. Creates a
    # temporary json file for the R script, and deletes it after use
    object, tf = tempfile.mkstemp()
    with open(tf, 'w') as f:
        json.dump(dmat, f, indent=4, sort_keys=True)
    distargs = ('Rscript', 'alleledistance.R',
                '--jsonfile', tf,
                '--outfile', os.path.join(outpath, outfile))
    subprocess.call(distargs)
    os.remove(tf)


def sequencetyping(dmat):
    # checks through the dictionary of different allele calls, those that have a score of 0 anywhere (identical
    # alleles) are assigned to the non-unique set, while those without a 0 score are assigned to the
    # unique set. Returns a set each of all unique and non-unique strains, and a dictionary of which strains
    # the non-unique strains match with (grouped together those that have a 0 score with each other).
    nonunique = set()
    unique = set()
    nu = collections.defaultdict(set)
    for strain1 in dmat:
        for strain2 in dmat[strain1]:
            if dmat[strain1][strain2] == 0:
                nonunique.add(strain1)
                nonunique.add(strain2)
                nu[strain1].add(strain2)
    for strain1 in dmat:
        if strain1 not in nonunique:
            unique.add(strain1)
    return unique, nonunique, nu


def stwriter(unique, nonunique, nu, outpath, stout):
    # writes out a csv file in the format of first row: summary information of number of strains and groups,
    # next rows: groups of strains
    # next rows: all remaining unique strains
    i = 0
    with open(os.path.join(outpath, stout), 'w') as f:
        file = csv.writer(f, delimiter = ' ')
        file.writerow(['number of unique strains: '+str(len(unique)),
                        'number of non-unique strains: '+str(len(nonunique)),
                        'number of non-unique groups: '+str(len(nu))])
        for strain1 in nu:
            row = ['group {}'.format(i), strain1]
            for strain2 in nu[strain1]:
                row.append(strain2)
            file.writerow(row)
            i += 1
        for x in unique:
            file.writerow(['unique', x])



def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', default='alleledistances.csv')
    parser.add_argument('--seqtyp', default='ST.csv')
    parser.add_argument('--outpath', default='mistreport/')
    parser.add_argument('path')
    return parser.parse_args()


def process(path, outpath, outfile, stout):
    pathfinder(outpath)
    files = filegetter(path)
    genedic = {}
    genedic2 = {}
    for file in files:
        reader(file, genedic, genedic2)
    dmat = compare(genedic, genedic2)
    unique, nonunique, nu = sequencetyping(dmat)
    csvwriter(dmat, outpath, outfile)
    stwriter(unique, nonunique, nu, outpath, stout)


def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.seqtyp)

if __name__ == '__main__':
    main()