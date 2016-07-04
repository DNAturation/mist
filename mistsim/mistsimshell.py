'''Main program that executes the simulation. Simulates effects of removing numbers of genes/genomes, and how many remain
Note: overwrites previously simulations within the same folder'''

import mistsim
import mistpopper
import mistchoppop
import subprocess
import argparse
import glob
import os
import multiprocessing

def grapher(outpath, outfile, outgraph):
    '''makes graph for the simulation of increasing threshold and how many genes/genomes are remaining'''
    garg = ('Rscript', 'mistgraph.R',
            outpath, outfile, outgraph)
    subprocess.call(garg)

def popgraph(outpath, outpop, popgraph):
    '''makes graph for simulation of chopping off genes and looking at how many genomes are perfect'''
    parg = ('Rscript', 'mistpopgraph.R',
            outpath, outpop, popgraph)
    subprocess.call(parg)


def clusters(outpath, chopout, outcut, startchop, outfile, corecalls):
    if not os.access(outcut, os.F_OK):
        os.mkdir(outcut)
    carg = ['Rscript', 'mistcluster.R', corecalls,
            '--chopsyms', os.path.join(outpath, chopout),
            '--chopped', outfile,
            '--outpath', outpath,
            '--outcuts', outcut,
            '--startchops'] + startchop
    subprocess.call(carg)



def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='path to report file')
    parser.add_argument('-o', '--outpath', default='/home/phac/kye/misty/sim/', help='directory for the summary file')
    subparsers = parser.add_subparsers(dest='subfunction')
    a_parser = subparsers.add_parser('Sim')
    a_parser.add_argument('--genethreshhold', default=9999999, help='maximum number within a gene of missing genomes tolerated before gene is removed')
    a_parser.add_argument('--genomethreshhold', default=9999999, help='maximum number within a genome of missing genes tolerated before genome is removed')

    a_parser.add_argument('--startpop', default=None, help='start cutoff of genes missing this many genomes')
    a_parser.add_argument('--endpop', default=None, help='end cutoff of genes missing this many genomes')
    a_parser.add_argument('--genomemin', default=0, type=int, help='minimum number cutoff to start at for genome')
    a_parser.add_argument('--genemin', default=0, type=int, help='minimum number cutoff to start at for gene')

    a_parser.add_argument('--outfile', default='simulator.csv', help='name of file to be created')
    a_parser.add_argument('--outpop', default='popper.csv', help='name of file to be created for popper')
    a_parser.add_argument('--outgraph', default='graph.png', help='name of output graph file')
    a_parser.add_argument('--poppedgraph', default='popper.png', help='name of output popper graph file')

    b_parser = subparsers.add_parser('Chop')
    b_parser.add_argument('--startpop', default=None, help='start cutoff of genes missing this many genomes')
    b_parser.add_argument('--endpop', default=None, help='end cutoff of genes missing this many genomes')
    b_parser.add_argument('--startchop', default=None, nargs='+')
    b_parser.add_argument('--outfile', default='chopped.csv')
    b_parser.add_argument('--outcut', default='cuts/')
    b_parser.add_argument('--chopout', default='chopsyms/')
    b_parser.add_argument('--alleles', default='/home/phac/kye/autocreate/prissy/alleles/')

    b_parser.add_argument('--distout', default='distances/')
    b_parser.add_argument('--corecalls', default='/home/phac/kye/autocreate/prissy/core_calls.csv')
    b_parser.add_argument('--outtree', default='distancegraph')

    b_parser.add_argument('--cores', default=multiprocessing.cpu_count())
    return parser.parse_args()


def processSim(path, outpath, outfile, genemin, genethreshhold, genomemin, genomethreshhold,
               outgraph, outpop, startpop, endpop, poppedgraph):
    print('Simulating threshold cutoffs...')
    mistsim.process(path, outpath, outfile, genemin, genethreshhold, genomemin, genomethreshhold)
    print('Graphing threshold csv...')
    grapher(outpath, outfile, outgraph)
    print('Popping genes...')
    mistpopper.process(path, outpath, outpop, startpop, endpop)
    print('Graphing popped csv...')
    popgraph(outpath, outpop, poppedgraph)

def processChop(path, outpath, outfile, outcut, startpop, endpop, startchop, chopout, alleles, distout, corecalls, outtree, cores):
    print('Chopping at given points')
    mistchoppop.process(path, outpath, outfile, startpop, endpop, startchop, chopout, alleles, cores)
    print('creating distance matrices and clustering genes...')
    clusters(outpath, chopout, outcut, startchop, outfile, corecalls)

def main():
    args = arguments()
    if args.subfunction == 'Sim':
        processSim(args.path, args.outpath, args.outfile, args.genemin, args.genethreshhold,
                   args.genomemin, args.genomethreshhold, args.outgraph, args.outpop,
                   args.startpop, args.endpop, args.poppedgraph)
    if args.subfunction == 'Chop':
        processChop(args.path, args.outpath, args.outfile, args.outcut, args.startpop, args.endpop,
                    args.startchop, args.chopout, args.alleles, args.distout, args.corecalls, args.outtree, args.cores)

if __name__ == '__main__':
    main()