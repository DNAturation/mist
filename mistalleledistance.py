# takes in a directory of MIST output files (jsons) and creates a distance calculation of the alleles.
# If two strains have similar alleles, their numbers (intersection point on CSV file) are lower. Counts similarity
# by checking which allele MIST matched each read with; different alleles increase the score
#
# Note that this script assigns 'NA' values when a match for a gene could not be found, and these values
# contribute to the score (matching with other NA's).
#
# Also groups strains depending on their allele matches;
# if their distance (number of different allele matches) are below a defined threshold, they are considered to be
# part of the same strain type.
#
# Note that this also combines strain groups together if they have a common strain that is in their strain types,
# even if individually they may cross the threshold
#
# (ex. strain1 and strain3 have a distance of 5, while strain1 and strain2 have a distance of 3, and strain3
# and strain2 have a distance of 3, at a set threshold of 3, all these will be grouped together despite the distance
# of 5 because strain2 connects strain1 and strain3)

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
    files = glob.glob(os.path.join(path, "*.json"))
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


def distget(allele1, allele2, dist):
    #increases distance score on difference in allele match
    if allele1 != allele2:
        dist += 1
    return dist


def compare(genedic, genedic2):
    # dmat is a dictionary that will be turned into a matrix.
    # It holds key of strain1, value of a dictionary(key = strain2, value = score)
    dmat = collections.defaultdict(dict)
    for strain1 in genedic:
        genedic2.pop(strain1)  # removes the strain from the second dictionary to prevent double comparison or two way comparison (A->B, B->A)
        for strain2 in genedic2:
            dist = 0  # distance counter
            if strain1 != strain2:
                for gene in genedic[strain1]:  # checks through each gene in the strains. This assumes that all strains have the same genes being checked
                    dist = distget(genedic[strain1][gene], genedic2[strain2][gene], dist)  # increases dist count by 1 on mismatch
                    dmat[strain1][strain2] = dist
    return dmat


def csvwriter(dmat, outpath, outfile, heatmap):
    # outsources the prettyifying of the dictionary into a table to an R script, alleledistance.R. Creates a
    # temporary json file for the R script, and deletes it after use
    object, tf = tempfile.mkstemp()
    with open(tf, 'w') as f:
        json.dump(dmat, f, indent=4, sort_keys=True)
    distargs = ('Rscript', 'alleledistance.R',
                '--jsonfile', tf,
                '--outfile', os.path.join(outpath, outfile),
                '--heatmap', os.path.join(outpath, heatmap))
    subprocess.call(distargs)
    os.remove(tf)


def grouping(nu, strain1, strain2, i):
    # groups non-unique strains together based on
    try:
        if strain1 in nu[i] or strain2 in nu[i]:
            nu[i].append(strain1)
            nu[i].append(strain2)
        else:
            i+=1
            grouping(nu, strain1, strain2, i)
    except KeyError:
        nu[i] = [strain1, strain2]


def combinegroup(nu):
    # takes in the groupings and makes sure that every group is actually different.
    # compares the strains within the groups, if there is the same strain name within two different groups,
    # the groups are joined together into one larger group
    iterd1 = {k: v for k, v in nu.items()}  # makes two dictionaries to loop over and compare with
    iterd2 = {k: v for k, v in nu.items()}
    newnu = {}
    used = []
    i=1
    change = False  # set flag to false for next iteration. Will be set to true on combining a group
    for loop1 in iterd1: # loops through the keys of the first dictionary
        iterd2.pop(loop1)  # removes the key form the second dictionary to prevent comparing the same key value pairs
        L1=True  # flag to decide whether or not to append the keys of the first loop to the new dictionary
                # appends if it wasn't already done so through combining groups
        for loop2 in iterd2:
            if len(iterd1[loop1].intersection(iterd2[loop2])) > 0:  # evaluates to true if any strain within the values (a set), matches the second values
                try:
                    newnu[i].update(iterd2[loop2])
                except KeyError:
                    newnu[i] = iterd1[loop1] | iterd2[loop2]  # combines the sets
                change=True  # flag set to show that a group was combined (= recurse through function)
                L1=False  # flag set to show that the key in loop1 does not need to be added again
            elif iterd2[loop2] not in used:
                newnu[i] = iterd2[loop2]  # if not a match, adds the set in alone without combining into a large group
                used.append(iterd2[loop2])
                i+=1
        if (L1) and (iterd1[loop1] not in used):  # if the first key's values were not added to the new dictionary(i.e. no matches), adds it in
            newnu[i] = iterd1[loop1]
            used.append(iterd1[loop1])
            i+=1
    return newnu, change


def sequencetyping(dmat, thresh):
    # checks through the dictionary of different allele calls, those that have a score of 0 anywhere (identical
    # alleles) are assigned to the non-unique set, while those without a 0 score are assigned to the
    # unique set. Returns a set each of all unique and non-unique strains, and a dictionary of which strains
    # the non-unique strains match with (grouped together those that have a 0 score with each other).
    nonunique = set()
    unique = set()
    nu = {}
    for strain1 in sorted(dmat, key=lambda k: len(dmat[k]), reverse=True):  # sorts from the longest distance matrix first to maximize chances of getting better matches initially
        for strain2 in dmat[strain1]:
            if dmat[strain1][strain2] <= int(thresh):  # if the allele distance is less than a defined amount, consider
                                                        # them to be part of the same strain
                nonunique.add(strain1)
                nonunique.add(strain2)  # adds both strains that have a perfect match to the non-unique group
                i=1
                grouping(nu, strain1, strain2, i)
    for strain1 in dmat:  # runs after nonunique set is completed
        if strain1 not in nonunique:  # if a strain isn't in the non-unique set, it must be unique
            unique.add(strain1)
    for group in nu:
        nu[group] = set(nu[group])
    change = True  # flag used for the combinegroup function on whether or not to rerun the function
    while change:
        nu, change = combinegroup(nu)  # calls function to combine all groups with the same strains until change == false

    return unique, nonunique, nu


def stwriter(unique, nonunique, nu, outpath, stout):
    # writes out a csv file in the format of first row: summary information of number of strains and groups,
    # next rows: groups of strains
    # next rows: all remaining unique strains
    with open(os.path.join(outpath, stout), 'w') as f:
        file = csv.writer(f, delimiter = ' ')
        file.writerow(['number of unique strains: '+str(len(unique)),
                        'number of non-unique strains: '+str(len(nonunique)),
                        'number of non-unique groups: '+str(len(nu))])  # summary information written to the top of file
        for group in nu:
            row = ['group {}:'.format(group)]  # counting of group numbers
            for strains in set(nu[group]):  # adds strains that are part of the same sequence type; belong to the same group
                row.append(strains)
            file.writerow(row)
        for x in unique:
            file.writerow(['unique', x])  # write each unique strain on its own row



def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', default='alleledistances.csv', help='name of output distance report csv file')
    parser.add_argument('--heatmap', default='distanceheatmap.png', help='name of the output file for the heatmap of alleledistances')
    parser.add_argument('--seqtyp', default='ST.csv', help='name of output sequence type report csv file')
    parser.add_argument('--outpath', default='mistreport/', help='output directory')
    parser.add_argument('--thresh', default=0, help='how many allele call differences are tolerated before assigning to different strains')
    parser.add_argument('path')
    return parser.parse_args()


def process(path, outpath, outfile, heatmap, stout, thresh):
    pathfinder(outpath)
    files = filegetter(path)
    genedic = {}
    genedic2 = {}
    for file in files:
        reader(file, genedic, genedic2)
    dmat = compare(genedic, genedic2)
    unique, nonunique, nu = sequencetyping(dmat, thresh)
    csvwriter(dmat, outpath, outfile, heatmap)
    stwriter(unique, nonunique, nu, outpath, stout)


def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.heatmap, args.seqtyp, args.thresh)

if __name__ == '__main__':
    main()