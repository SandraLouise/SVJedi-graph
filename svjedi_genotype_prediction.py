#!/usr/bin/env python3

"""*******************************************************************************
    Adapted from SVJedi/modules/genotype.py v.1.1.4 to work with json input.
    Git: https://github.com/llecompte/SVJedi
*******************************************************************************"""

import sys
import numpy as np
import argparse
import math
from decimal import *
import json

def main(args):
    """ Parsing arguments """
    parser = argparse.ArgumentParser(description="Structural variations genotyping using long reads")
    
    parser.add_argument("-d", "--aln", metavar="<alndict>", nargs=1, required=True)
    
    parser.add_argument("-v", "--vcf", metavar="<vcffile>", help="vcf format", required=True)
    
    parser.add_argument("-o", "--output", metavar="<output>", nargs=1, help="output file")

    parser.add_argument("-ladj", metavar="<allele_size>", nargs=1, type=int, default=[5000], help="Sequence allele adjacencies at each side of the SV")

    parser.add_argument(
        "-ms",
        "--minsupport",
        metavar="<minNbAln>",
        type=int,
        default=3,
        help="Minimum number of alignments to genotype a SV (default: 3>=)",
    )

    args = parser.parse_args()

    # check if outputfile
    if args.output is None:
        output = "genotype_results.txt"
    else:
        output = args.output[0]

    vcffile = args.vcf
    aln_dict_file = args.aln[0]
    min_support = args.minsupport
    l_adj = args.ladj[0]

    with open(aln_dict_file, 'r') as file:
        dict_of_informative_aln = json.load(file)


    decision_vcf(dict_of_informative_aln, vcffile, output, min_support, l_adj)

def decision_vcf(dictReadAtJunction, inputVCF, outputDecision, minNbAln, l_adj):
    """ Output in VCF format and take genotype decision """
    getcontext().prec = 28
    outDecision = open(outputDecision, "w")
    with open(inputVCF) as inputFile:
        for line in inputFile:
            if line.startswith("##"):
                outDecision.write(line)

            elif line.startswith("#C"):
                outDecision.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                outDecision.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">\n')
                outDecision.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Depth of each allele by sample">\n')
                outDecision.write(line.rstrip("\n") + "\t" + "\t".join(["FORMAT", "SAMPLE"]) + "\n")

            else:
                in_chrom, in_start, _, __, in_type, ___, ____, in_info, *_ = line.rstrip("\n").split("\t")

                ### get SVTYPE ###
                if 'SVTYPE' in in_info:
                    if in_info.split(';')[-1].startswith('SVTYPE='):
                        svtype = in_info.split('SVTYPE=')[1]
                    else:
                        svtype = in_info.split('SVTYPE=')[1].split(';')[0]
                else:
                    svtype = ''
                
                ### get LENGTH for DELETION ###                 
                if svtype == 'DEL':
                    if "SVLEN=FALSE" in in_info:
                        if in_info.startswith("END="):
                            end = in_info.split("END=")[1].split(";")[0]
                        else:
                            end = in_info.split(";END=")[1].split(";")[0]

                        in_length = int(end) - int(in_start)

                    elif "SVLEN=" in in_info:
                        if in_info.split(';')[-1].startswith('SVLEN='):
                            in_length = abs(int(in_info.split("SVLEN=")[1]))
                        else:
                            in_length = abs(int(in_info.split("SVLEN=")[1].split(";")[0]))
                
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV
                    
                
                ### get LENGTH for INSERTION ###                    
                elif svtype == 'INS': 
                    in_length = len(in_type) 
                    
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV
                
                
                ### get LENGTH for INVERSION ###                
                elif svtype == 'INV':
                    if in_info.startswith("END="):
                        end = in_info.split("END=")[1].split(';')[0]
                    else:
                        end = in_info.split(";END=")[1].split(';')[0]               
                    in_length = int(end) - int(in_start)
                    
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV


                ### get sv id for TRANSLOCATION ###             
                elif svtype == 'BND': 
                    if in_info.startswith("END="):
                        end = in_info.split("END=")[1].split(';')[0]
                    elif ";END=" in in_info:
                        end = in_info.split(";END=")[1].split(';')[0]   
                    elif "[" in  in_type:
                        end = in_type.split(':')[1].split('[')[0]
                    elif "]" in  in_type:
                        end = in_type.split(':')[1].split(']')[0]
                    
                    if 'CHR2=' in in_info: 
                        chr2 = in_info.split('CHR2=')[1].split(';')[0]
                    elif '[' in in_type:
                            chr2 = in_type.split(':')[0].split('[')[1]     #ALT[CHR2:POSTION[
                    elif ']' in in_type:
                            chr2 = in_type.split(':')[0].split(']')[1]     #ALT]CHR2:POSTION]
                            
                    in_sv = in_chrom + "_" + in_start + "-" + chr2 + "-" + end #define sv id for TRANS
                
                
                #######################################################################################
                #Asign genotype 
                
                if svtype in ('DEL', 'INS', 'INV', 'BND') and in_sv in list(dictReadAtJunction.keys()) and abs(in_length) >= 50:
                    nbAln = [len(x) for x in dictReadAtJunction[in_sv]]
                    genotype, proba = likelihood(nbAln, svtype, in_length, minNbAln, l_adj)
            
                else: #if svtype different from DEL, INS, INV, BND or if sv not in supported by alignment
                    nbAln = [0,0]
                    genotype = "./."
                    proba = [".",".","."]
                
                    
                #######################################################################################
                #Output genotype in VCF
                
                numbers = ",".join(str(y) for y in nbAln)
                if len(line.split("\t")) <= 8:
                    new_line = (
                        line.rstrip("\n")
                        + "\t"
                        + "GT:DP:AD:PL"
                        + "\t"
                        + genotype
                        + ":"
                        + str(round(sum(nbAln), 3))
                        + ":"
                        + str(numbers)
                        + ":"
                        + str(','.join(proba))
                    )
                    outDecision.write(new_line + "\n")

                else:
                    line_without_genotype = line.split("\t")[0:8]
                    new_line = (
                        "\t".join(line_without_genotype)
                        + "\t"
                        + "GT:DP:AD:PL"
                        + "\t"
                        + genotype
                        + ":"
                        + str(round(sum(nbAln), 3))
                        + ":"
                        + str(numbers)
                        + ":"
                        + str(','.join(proba))
                    )
                    outDecision.write(new_line + "\n")

    outDecision.close()

def likelihood(all_count, svtype, svlength, minNbAln, l_adj):
    """ Compute likelihood """
    
    unbalanced_sv = ("DEL", "INS")
    if svtype in unbalanced_sv:
        c1, c2 = allele_normalization(all_count, svtype, svlength, l_adj)  # allelic sequence normalization for unbalanced SV
    else:
        c1, c2 = all_count
    
    rc1 = int(round(c1,0))
    rc2 = int(round(c2,0))
    e = 0.00005 #sequencing err

    lik0 = Decimal(c1*math.log10(1-e)) + Decimal(c2*math.log10(e)) 
    lik1 = Decimal((c1+c2)*math.log10(1/2)) 
    lik2 = Decimal(c2*math.log10(1-e)) + Decimal(c1*math.log10(e))
    
    L = [lik0, lik1, lik2]
    
    index_of_L_max = [i for i, x in enumerate(L) if x == max(L)]
    if len(index_of_L_max) == 1: 
        geno_not_encoded = str(index_of_L_max[0])
        geno = encode_genotype(geno_not_encoded)
        
    else:
        geno = "./."    #no genotype estimation since likelihood are not conclusive
    
    #Check for minimum number of alignment to assign a genotype
    if not sum(all_count) >= minNbAln: # minimum support
        geno = "./."
                
    combination = Decimal(math.log10(math.comb(rc1 + rc2, rc1)))
    lik0 += combination 
    lik1 += combination
    lik2 += combination

    #phred scaled score
    prob0 = -10*lik0
    prob1 = -10*lik1
    prob2 = -10*lik2                                 
                 
    prob = [str(int(prob0)), str(int(prob1)), str(int(prob2))]
            
    return geno, prob

def allele_normalization(nb_aln_per_allele, svtype, svlength, l_adj):
    ''' Allele length normalization '''
    if svlength > (2*l_adj): svlength = 2*l_adj #for upper bound, if case of sv size > 2XLadj, only 2 sequences of 2Ladj are represented
    
    if svtype == "DEL":
        nb_aln_longest_allele_seq = nb_aln_per_allele[0]
        if nb_aln_longest_allele_seq > 0:
            nb_aln_per_allele[0] = round(nb_aln_longest_allele_seq * (2*l_adj) / ((2*l_adj) + svlength), 3)
    
    elif svtype == "INS":
        nb_aln_longest_allele_seq = nb_aln_per_allele[1]
        if nb_aln_longest_allele_seq > 0:
            nb_aln_per_allele[1] = round(nb_aln_longest_allele_seq * (2*l_adj) / ((2*l_adj) + svlength), 3)
    
    return nb_aln_per_allele

def encode_genotype(g):
    ''' Encode genotype from 0, 1, 2 to 0/0, 0/1, 1/1 '''
    if g == '0': genotype = "0/0"
    elif g == '1': genotype = "0/1"
    elif g == '2': genotype = "1/1"
    
    return genotype

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])