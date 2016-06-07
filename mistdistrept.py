#takes input from parsed JSON files from MIST and sets threshhold
#report part1 = dictionary of strains, values of genes they are missing
#report part2 = dictionary of genes, values of strains that lack those genes

import argparse
import json
import glob
import os
import re

def pathfinder(outpath):
    if not os.access(outpath, os.F_OK):
        os.mkdir(outpath)

def fileget(path):
    files = glob.glob(path+'*.json')
    return files

def reader(files):
    with open (files, 'r') as f:
        data = json.load(f)
    return data

# def repper(outpath, outfile):
#     with open(os.path.join(outpath, outfile), 'r') as f:
#         rep = f.read()
#     return rep

def writer(data, testtype, files):
    '''previously wrote data to files (hence the name), now provides the strain name,
    number of genes the genome does not have, and a list of genes that did not match for the genome
    used as primary source of information for part1 of report'''
    genes = data['Results'][0]['TestResults'][testtype]
    genesmissing = 0
    genelist = []
    for gene in genes:
        # print (genes[gene]['CorrectMarkerMatch'] == False)
        if genes[gene]['CorrectMarkerMatch'] == False:
            genelist.append(gene)
            genesmissing+=1
    # dmissinggenes[os.path.splitext(os.path.basename(files))[0][:-len(testtype)]] = \
    #     (str(genelist).strip('[]') + ' total = ' + str(genesmissing))
    return os.path.splitext(os.path.basename(files))[0][:-len(testtype)], genesmissing, genelist

#count number of genomes a gene is missing from
# def misscounter(data, testtype, rep):
#     genes = data['Results'][0]['TestResults'][testtype]
#     mostmisses = 0
#     missinggene = []
#     dmissingfrom = {}
#     for gene in genes:
#         misses = re.findall('{}'.format(gene), str(rep).strip('[]'))
#         # f.write(gene+' missing from '+str(len(misses))+' genomes \n')
#         if len(misses) > mostmisses:
#             mostmisses = len(misses)
#             missinggene = [gene]
#         elif len(misses) == mostmisses:
#             missinggene.append(gene)
#     for x in missinggene:
#         dmissingfrom[x]=mostmisses
#     return dmissingfrom


# def theouts(d):
#     genomemisses = 0
#     mostmissinggenomes = []
#     dgenomes = {}
#     for genome in d:
#         if d[genome] > genomemisses:
#             genomemisses = d[genome]
#             mostmissinggenomes = [genome]
#         elif d[genome] == genomemisses:
#             mostmissinggenomes.append (genome)
#     for x in mostmissinggenomes:
#         dgenomes[x]=genomemisses
#     return dgenomes

def genes(genelist, dwriter):
    '''takes in a list of all genes that have misses and a dictionary of k=strain v=list of missing genes,
    returns a dictionary of k=genes v=strains that lack it for part 2 of the report'''
    dmisslist={}
    for gene in genelist:
        dmisslist[gene]=[]
    for strain in dwriter:
        for index, genes in enumerate(dwriter[strain]):
            dmisslist[genes].append(strain)
    return dmisslist



def JSONwriter(outpath, outfile, dwriter, dmisslist):
    '''writes the report into a json file'''
    with open(os.path.join(outpath, outfile), 'a') as f:
        data = {'GenomesMissingGenes': dwriter, 'GenesMissingGenomes': dmisslist}
        json.dump(data, f, indent=4, sort_keys=True)

def genetotal(mistfolder, testtype):
    with open(os.path.join(mistfolder, testtype)+'.markers', 'r') as f:
        data = f.readlines()
        genelist = []
        for x in range(1, len(data)):
            genelist.append(data[x].split()[0])
        return genelist

        # genelist = re.findall('\w*\.fasta', data)
        # for gene in genelist:
        #     os.path.splitext(gene)
        # return genelist



def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outpath', default='./mistfail/')
    parser.add_argument('--outfile', default='proteinsfailed.json')
    parser.add_argument('-t', '--testtype', default='MLST')
    parser.add_argument('--mistfolder', default='/home/cintiq/Desktop/campylobacterjejuni/')
    parser.add_argument('path')
    return parser.parse_args()

def process(path, outpath, outfile, testtype, mistfolder):
    pathfinder(outpath)
    files = fileget(path)
    dwriter = {}
    dgenomes={}
    dmostmissing={}
    d = {}
    dmisslist={}
    # rep=[]
    genelist = genetotal(mistfolder, testtype)
    for item in files:
        # dgenometotal={}
        data = reader(item)
        missingno = writer(data, testtype, item)
        # rep.extend(missingno[2])
        d[missingno[0]]=missingno[1]
        # dgenometotal['Missing'] = missingno[2]
        # # print (dgenometotal)
        # dwriter[item]=dgenometotal
        dwriter[missingno[0]] = missingno[2]
        #dwriter = of each genome (file), what genes are missing
    # print(dgenometotal)
    # dmisslist['GenesMissingFromGenome'] = dwriter
    # dmostmissing['MostMissingGenes'] = misscounter(reader(files[0]), testtype, rep)
    #dmostmissing = which genes have highest occurrence of missing from genomes
    # dgenomes['LeastMatchingGenomes'] = theouts(d)
    #dgenomes = genomes that have the least matches
    dmisslist=genes(genelist, dwriter)
    JSONwriter(outpath, outfile, dwriter, dmisslist)

def main():
    args = arguments()
    process(args.path, args.outpath, args.outfile, args.testtype, args.mistfolder)


if __name__ == '__main__':
    main()