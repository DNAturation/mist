#takes input from parsed JSON files from MIST and sets threshhold, copies files that pass

import argparse
import json
import glob
import shutil
import os

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

def counter(data, testtypename):
    genes = data['Results'][0]['TestResults'][testtypename]
    missingno = 0
    for gene in genes:
        if genes[gene]['CorrectMarkerMatch'] == False:
            missingno+=1
        else:
            continue

    return missingno


def cull(file, outpath, missing, threshold):
    if missing <= threshold:
        strain = os.path.basename(file)
        shutil.copy(file, outpath + strain)
    else:
        pass




def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-thresh', '--threshhold', type=int, required=True)
    parser.add_argument('-o', '--outpath', default='./mistpass/')
    parser.add_argument('--testtypename', default='core')
    parser.add_argument('path')
    return parser.parse_args()

def process(path, outpath, testtypename, threshhold):
    pathfinder(outpath)
    files = fileget(path)
    for item in files:
        data = reader(item)
        missing = counter(data, testtype)
        cull(item, outpath, missing, threshhold)

def main():
    args = arguments()
    process(args.path, args.outpath, args.testtypename, args.threshhold)


if __name__ == '__main__':
    main()