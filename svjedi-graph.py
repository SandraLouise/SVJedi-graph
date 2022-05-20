#!/usr/bin/env python3

"""*******************************************************************************
    Name: SVJedi-graph
    Description: SVjedi-graph aims to genotype structural variant with long reads data using a variation graph.
    Author: Sandra Romain
    Contact: sandra.romain@inria.fr, Inria/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
    
    Copyright (C) 2019 Inria
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""

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
    threads = args.threads

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
