# takes in a directory of MIST output files (jsons) and creates a distance calculation of the alleles.
# If two strains have similar alleles, their numbers (intersection point on CSV file) are lower. Counts similarity
# by checking which allele MIST matched each read with; different alleles increase the score
import json
import argparse
import os
import glob
import subprocess
import collections
import tempfile


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


def csvwriter(dmat, outfile):
    # outsources the prettyifying of the dictionary into a table to an R script, alleledistance.R. Creates a
    # temporary json file for the R script, and deletes it after use
    tf = tempfile.NamedTemporaryFile()
    with open(tf, 'w') as f:
        json.dump(dmat, f, indent=4, sort_keys=True)
    distargs = ('Rscript', 'alleledistance.R',
                '--jsonfile', tf,
                '--outfile', outfile)
    subprocess.call(distargs)
    os.remove(tf)



def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', default='alleledistances.csv')
    parser.add_argument('path')
    return parser.parse_args()


def process(path, outfile):
    files = filegetter(path)
    genedic = {}
    genedic2 = {}
    for file in files:
        reader(file, genedic, genedic2)
    dmat = compare(genedic, genedic2)
    csvwriter(dmat, outfile)


def main():
    args = arguments()
    process(args.path, args.outfile)

if __name__ == '__main__':
    main()