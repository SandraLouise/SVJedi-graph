#!/usr/bin/env python3

import sys
import argparse

def main(args):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v", 
        "--vcf", 
        metavar="<inputVCF>", 
        type=str,
        nargs=1, 
        required=True)

    parser.add_argument(
        "-r", 
        "--ref", 
        metavar="<referenceGenome>", 
        type=str,
        nargs=1, 
        required=True)

    parser.add_argument(
        "-l", 
        "--ladj", 
        metavar="<lengthOfFlankingRegions", 
        type=int,
        nargs=1, 
        default=5000)
    
    parser.add_argument(
        "-n", 
        "--number", 
        metavar="<numberOfExtractedSV", 
        type=int,
        nargs=1, 
        required=False)
    
    parser.add_argument(
        "-c", 
        "--cut", 
        action='store_true')
    
    parser.add_argument(
        "-e", 
        "--expid", 
        metavar="<experienceID", 
        type=str,
        required=False)

    args = parser.parse_args()

    inputVCF = args.vcf[0]
    reference_fasta = args.ref[0]
    ladj = int(args.ladj)

    if args.cut:
        cut = True
    else:
        cut = False
    
    ########### Tests
    if args.number:
        n = int(args.number[0])
    
    if args.expid:
        newRef_file_name = args.expid + '_' + 'reference_subsequences.fa'
        newVCF_file_name = args.expid + '_' + 'converted_variants.vcf'
    ###########
    
    else:
        newRef_file_name = 'reference_subsequences.fa'
        newVCF_file_name = 'converted_variants.vcf'


    #1. Get reference genome sequence

    dict_of_chrom = {}
    with open(reference_fasta) as sequenceFile:
        sequence = ""
        for line in sequenceFile:
            if line.startswith(">"):
                if sequence != "":
                    dict_of_chrom[header] = sequence
                sequence = ""
                header = line.rstrip("\n")[1:].split()[0]

            else:
                sequence += line.rstrip("\n")

        dict_of_chrom[header] = sequence


    #2. Read SV data

    count = 0

    with open(inputVCF) as vcffile:
        vcf_header = ""
        newRef_file = open(newRef_file_name, "w")
        newVCF_file = open(newVCF_file_name, "w")

        for line in vcffile:

            if line[0] == "#":
                vcf_header = vcf_header + line

            else:
                
                #Test
                if args.number:
                    if count == n:
                        break

                contig, pos, sv_id, __, __, __, __, info, __, __ = line.rstrip().split("\t")

                #Get SV type
                if 'SVTYPE' in info:
                    if info.split(';')[-1].startswith('SVTYPE='):
                        svtype = info.split('SVTYPE=')[1]
                    else:
                        svtype = info.split('SVTYPE=')[1].split(';')[0]
                # else:
                #     svtype = ''

                #Get breakpoint positions on ref
                start = int(pos)

                if info.split(';')[-1].startswith('SVLEN='):
                    length = abs(int(info.split("SVLEN=")[1]))
                else:
                    length = abs(int(info.split("SVLEN=")[1].split(";")[0]))


                if svtype == 'DEL':
                    end = start + length 
                
                # elif svtype == 'INT':
                #     end = start #ou start + 1 ?

                #Get ref seq
                new_contig_id = "_".join(["ref", contig, "-".join([pos, str(length)])])
                newRef_file.write(">" + new_contig_id + "\n")

                if cut and (length > 2 * ladj):
                    left_part = dict_of_chrom[contig][int(pos) - ladj + 1 : int(pos) + ladj + 1]
                    # print(len(left_part))
                    right_part = dict_of_chrom[contig][int(end) - (ladj - 1) : int(end) + ladj + 1]
                    # print(len(right_part))
                    newRef_file.write(left_part + right_part + "\n")

                else:
                    newRef_file.write(dict_of_chrom[contig][int(pos) - ladj + 1 : int(end) + ladj + 1] + "\n")

                #Create new VCF with updated positions
                if count == 0:
                    newVCF_file.write(vcf_header)

                newVCF_file.write(convert_sv(line, ladj, new_contig_id, cut))

                count += 1
                
        
        newRef_file.close()
        newVCF_file.close()


def convert_sv(sv_line, ladj, new_contig_id, cut):
    contig, pos, sv_id, ref, alt, qual, sv_filter, info, gt_format, truth = sv_line.split("\t")

    if info.split(';')[-1].startswith('SVLEN='):
        length = abs(int(info.split("SVLEN=")[1]))
    else:
        length = abs(int(info.split("SVLEN=")[1].split(";")[0]))

    new_pos = ladj

    if (length > 2 * ladj) and cut:
        new_end = 3 * ladj

        #Replace SVLEN value
        original_svlen_with_flag = 'SVLEN=' + '-' + str(length) #only for DEL
        new_svlen_with_flag = 'SVLEN=' + '-' + str(2 * ladj) #only for DEL
        info = info.replace(original_svlen_with_flag, new_svlen_with_flag)

    else:
        new_end = ladj + length

    #Replace END value
    if info.split(';')[0].startswith('END='):
        #other info flags can contain "END=" ("CIEND=")
        original_end = abs(int(sv_line.split("\tEND=")[1].split(";")[0]))
        original_end_with_flag = "\tEND=" + str(original_end) + ";"
        new_end_with_flag = "\tEND=" + str(new_end) + ";"
        
        new_sv_line = new_contig_id + "\t" + str(new_pos) + "\t" + sv_id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + sv_filter + "\t" + info + "\t" + gt_format + "\t" + truth
        new_sv_line.replace(original_end_with_flag, new_end_with_flag)
        return new_sv_line

    elif info.split(';')[-1].startswith('END='):
        original_end = abs(int(info.split("END=")[1]))
        original_end_with_flag = ";END=" + str(original_end)
        new_end_with_flag = ";END=" + str(new_end)

    else:
        original_end = abs(int(info.split(";END=")[1].split(";")[0]))
        original_end_with_flag = ";END=" + str(original_end) + ";"
        new_end_with_flag = ";END=" + str(new_end) + ";"

    new_info = info.replace(original_end_with_flag, new_end_with_flag)
    new_sv_line = new_contig_id + "\t" + str(new_pos) + "\t" + sv_id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + sv_filter + "\t" + new_info + "\t" + gt_format + "\t" + truth

    return new_sv_line

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])