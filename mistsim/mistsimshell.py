'''Main program that executes the simulation. Simulates effects of removing numbers of genes/genomes, and how many remain
Note: overwrites previously simulations within the same folder'''

import mistsim
import mistpopper
import subprocess
import argparse

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


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genethreshhold', default=9999999, help='maximum number within a gene of missing genomes tolerated before gene is removed')
    parser.add_argument('--genomethreshhold', default=9999999, help='maximum number within a genome of missing genes tolerated before genome is removed')

    parser.add_argument('--startpop', default=None, help='start cutoff of genes missing this many genomes')
    parser.add_argument('--endpop', default=None, help='end cutoff of genes missing this many genomes')
    parser.add_argument('--genomemin', default=0, type=int, help='minimum number cutoff to start at for genome')
    parser.add_argument('--genemin', default=0, type=int, help='minimum number cutoff to start at for gene')

    parser.add_argument('--outfile', default='simulator.csv', help='name of file to be created')
    parser.add_argument('--outpop', default='popper.csv', help='name of file to be created for popper')
    parser.add_argument('-o', '--outpath', default='./sim/', help='directory for the summary file')
    parser.add_argument('--outgraph', default='graph.png', help='name of output graph file')
    parser.add_argument('--popgraph', default='popper.png', help='name of output popper graph file')
    parser.add_argument('path', help='path to report file')
    return parser.parse_args()


def main():
    args = arguments()
    print('Simulating threshold cutoffs...')
    mistsim.process(args.path, args.outpath, args.outfile, args.genemin, args.genethreshhold, args.genomemin, args.genomethreshhold)
    print('Graphing threshold csv...')
    grapher(args.outpath, args.outfile, args.outgraph)
    print('Popping genes...')
    mistpopper.process(args.path, args.outpath, args.outpop, args.startpop, args.endpop)
    print('Graphing popped csv...')
    popgraph(args.outpath, args.outpop, args.popgraph)



if __name__ == '__main__':
    main()