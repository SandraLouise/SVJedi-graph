#!/usr/bin/env python3

import sys
import argparse
import subprocess

def main(args):

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", 
        "--vcf", 
        metavar="<inputVCF>", 
        type=str,
        required=True)
    parser.add_argument(
        "-r", 
        "--ref", 
        metavar="<referenceGenome>", 
        type=str, 
        required=True)
    parser.add_argument(
        "-q", 
        "--reads", 
        metavar="<queryReads>", 
        type=str, 
        required=True)
    parser.add_argument(
        "-p", 
        "--prefix", 
        metavar="<outFilesPrefix>", 
        type=str, 
        required=True)
    parser.add_argument(
        "-t", 
        "--threads", 
        metavar="<threadNumber>", 
        type=int, 
        default=[1])

    args = parser.parse_args()
    inVCF = args.vcf
    inREF = args.ref
    inFQ = args.reads
    outPrefix = args.prefix
    threads = args.threads[0]

    ladj = 5000

    #### Create variant graph
    outGFA = outPrefix + ".gfa"
    c1 = "python3 construct-graph.py -v {} -r {} -o {}".format(inVCF, inREF, outGFA)
    subprocess.run(c1, shell=True)

    #### Map reads on graph
    outGAF = outPrefix + ".gaf"
    aln_log = outGAF + ".log"
    c2 = "GraphAligner -g {} -f {} -a {} -x vg -t {} > {}".format(outGFA, inFQ, outGAF, threads, aln_log)
    subprocess.run(c2, shell=True)

    #### Filter alns
    outJSON = outPrefix + "_informative_aln.json"
    c3 = "python3 filter-alignments.py -a {} -g {} -p {}".format(outGAF, outGFA, outPrefix)
    subprocess.run(c3, shell=True)

    #### Genotype
    outVCF = outPrefix + "_genotype.vcf"
    c4 = "python3 predict-genotype.py -d {} -v {} -o {}".format(outJSON, inVCF, outVCF)
    subprocess.run(c4, shell=True)

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])