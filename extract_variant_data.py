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
    
    # parser.add_argument(
    #     "-n", 
    #     "--number", 
    #     metavar="<numberOfExtractedSV", 
    #     type=int,
    #     nargs=1, 
    #     required=False)
    
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

    parser.add_argument(
        "-V", 
        "--version", 
        metavar="<versionID", 
        type=str,
        required=False)

    args = parser.parse_args()

    inputVCF = args.vcf[0]
    reference_fasta = args.ref[0]
    l_adj = int(args.ladj)

    if args.cut:
        cut = True
    else:
        cut = False
    
    ########### Tests
    # if args.number:
    #     n = int(args.number[0])
    
    if args.expid and args.version:
        newRef_file_name = args.expid + '_' + 'reference_subsequences' + '_' + args.version + '.fa'
        newVCF_file_name = args.expid + '_' + 'converted_variants' + '_' + args.version + '.vcf'
    
    elif args.expid:
        newRef_file_name = args.expid + '_' + 'reference_subsequences.fa'
        newVCF_file_name = args.expid + '_' + 'converted_variants.vcf'
    ###########
    
    else:
        newRef_file_name = 'reference_subsequences.fa'
        newVCF_file_name = 'converted_variants.vcf'
    
    sv_min_dist = 5000


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


    #2. Find close SVs

    with open(inputVCF) as file: 

        prev_end = None
        dict_of_close_sv = dict()
        found_close_sv = False
        #dict = {chrom : {first_closeSV : [following close SVs]}}

        for line in file:
            if not line.startswith('#'):
                chrom, pos, sv_id, __, __, __, __, info, *__ = line.rstrip().split("\t")
                start = int(pos)
                end = start + get_sv_len(info) #only for deletions

                if chrom not in dict_of_close_sv:
                    dict_of_close_sv[chrom] = dict()

                if prev_end is not None:
                    if (start < prev_end + sv_min_dist) and (chrom == prev_chrom):
                        
                        #First close SV
                        if found_close_sv == False:
                            first_close_sv = prev_sv_id
                            dict_of_close_sv[chrom][first_close_sv] = [sv_id]
                            found_close_sv = True
                        
                        else:
                            dict_of_close_sv[chrom][first_close_sv].append(sv_id)
                    
                    #End of close SV 
                    else:
                        found_close_sv = False
                
                prev_sv_id = sv_id
                prev_chrom = chrom
                prev_end = end

    # print(dict_of_close_sv)

    #3. Extract and convert SV data

    with open(inputVCF) as file:
        newRef_file = open(newRef_file_name, "w")
        newVCF_file = open(newVCF_file_name, "w")
        found_close_sv = False
        prev_end = None

        for line in file:

            if line.startswith('#'):
                newVCF_file.write(line) 

            else:
                
                #Test
                # if args.number:
                #     if count == n:
                #         break

                #Get SV info
                chrom, pos, sv_id, __, __, __, __, info, *__ = line.rstrip().split("\t")
                start_on_chr = int(pos)
                length = get_sv_len(info)
                sv_type = get_sv_type(info)

                if sv_type == 'DEL':
                    end_on_chr = start_on_chr + length #only for deletions
                
                # elif sv_type == 'INS':

                elif sv_type == 'INV':
                    end_on_chr = int(info.split(';END=')[1].split(';')[0])

                #First SV of sequence of close SVs
                if sv_id in dict_of_close_sv[chrom]:
                    found_close_sv = True
                    first_close_sv = sv_id

                    #Get left flanking region + sv ref seq
                    region_id = '_'.join(['ref', chrom, '-'.join([sv_type, pos, str(length)])]) #= 'ref_' + chrom + '_' + sv_type + "-" + pos_sv_1 + '-' + len_sv_1 (ex: 'ref_1_DEL-88749-1400')
                    left_flanking_region = get_sv_reference(dict_of_chrom, chrom, start_on_chr - l_adj, start_on_chr, l_adj, False)
                    sv_refSeq = get_sv_reference(dict_of_chrom, chrom, start_on_chr, end_on_chr, l_adj, cut)
                    extended_refSeq = left_flanking_region + sv_refSeq
                    start_of_extended_refSeq = start_on_chr - l_adj

                    #Convert SV coordinates
                    converted_start = l_adj
                    if sv_type == 'DEL':
                        converted_end = converted_start + length #DELETION
                    elif sv_type == 'INV':
                        converted_end = converted_start + (end_on_chr - start_on_chr)
                    converted_vcf_lines = convert_sv(line, 'region_id', converted_start, converted_end, l_adj, cut)

                #Following close SVs
                elif found_close_sv and (sv_id in dict_of_close_sv[chrom][first_close_sv]):
                    #Get left flanking region + sv ref seq
                    region_id = '_'.join([region_id, '-'.join([sv_type, pos, str(length)])]) #ex: 'ref_1_DEL-88749-1400' + '_' + sv2_type + '-' + pos_sv_2 + '-' + len_sv_2

                    if start_on_chr <= prev_end:
                        #no flanking region after last sv and before current sv
                        sv_refSeq = get_sv_reference(dict_of_chrom, chrom, prev_end, end_on_chr, l_adj, cut)
                        extended_refSeq = extended_refSeq + sv_refSeq

                    else:
                        left_flanking_region = get_sv_reference(dict_of_chrom, chrom, prev_end, start_on_chr, l_adj, False)
                        sv_refSeq = get_sv_reference(dict_of_chrom, chrom, start_on_chr, end_on_chr, l_adj, cut)
                        extended_refSeq = extended_refSeq + left_flanking_region + sv_refSeq

                    #Convert SV coordinates
                    converted_start = start_on_chr - start_of_extended_refSeq
                    if sv_type == 'DEL':
                        converted_end = converted_start + length #DELETION
                    elif sv_type == 'INV':
                        converted_end = converted_start + (end_on_chr - start_on_chr)
                    converted_vcf_lines = converted_vcf_lines + convert_sv(line, 'region_id', converted_start, converted_end, l_adj, cut)
                    
                    #Last close SV
                    if sv_id == dict_of_close_sv[chrom][first_close_sv][-1]:
                        #Get right flanking region and write ref seq and vcf lines
                        right_flanking_region = get_sv_reference(dict_of_chrom, chrom, prev_end, prev_end + l_adj, l_adj, False)
                        extended_refSeq = extended_refSeq + right_flanking_region
                        newRef_file.write('>' + region_id + "\n")
                        newRef_file.write(extended_refSeq + "\n")

                        converted_vcf_lines = converted_vcf_lines.replace('region_id', region_id)
                        newVCF_file.write(converted_vcf_lines)

                        found_close_sv = False

                #SV not in sequence of close SV
                else:
                    region_id = '_'.join(['ref', chrom, '-'.join([sv_type, pos, str(length)])])
                    converted_start = l_adj

                    if sv_type == 'DEL':
                        converted_end = converted_start + length #DELETION
                    elif sv_type == 'INV':
                        converted_end = converted_start + (end_on_chr - start_on_chr)

                    newVCF_file.write(convert_sv(line, region_id, converted_start, converted_end, l_adj, cut))
                    newRef_file.write('>' + region_id + "\n")
                    newRef_file.write(dict_of_chrom[chrom][start_on_chr - l_adj + 1 : end_on_chr + l_adj + 1] + "\n")

                prev_end = end_on_chr 

        newRef_file.close()
        newVCF_file.close()

def get_sv_len(sv_info):

    if sv_info.split(';')[-1].startswith('SVLEN='):
        return abs(int(sv_info.split("SVLEN=")[1]))
    else:
        return abs(int(sv_info.split("SVLEN=")[1].split(";")[0]))

def get_sv_type(sv_info):
    if 'SVTYPE' in sv_info:
        if sv_info.split(';')[-1].startswith('SVTYPE='):
            return sv_info.split('SVTYPE=')[1]
        else:
            return sv_info.split('SVTYPE=')[1].split(';')[0]

def get_sv_reference(dict_of_chrom, chrom, sv_start, sv_end, l_adj, cut):
    if (sv_end - sv_start > 2 * l_adj) and cut:
        left_ref = dict_of_chrom[chrom][sv_start + 1 : sv_start + l_adj + 1]
        right_ref = dict_of_chrom[chrom][sv_end - (l_adj - 1) : sv_end + 1]
        return left_ref + right_ref

    else:
        return dict_of_chrom[chrom][sv_start + 1 : sv_end + 1]

def convert_sv(sv_line, region_id, converted_start, converted_end, l_adj, cut):
    __, __, sv_id, ref, alt, qual, sv_filter, info, *__ = sv_line.split("\t")
    length = converted_end - converted_start

    if (length > 2 * l_adj) and cut:
        converted_end = 3 * l_adj
        #Replace SVLEN value
        original_svlen_with_flag = 'SVLEN=' + '-' + str(length) #only for DEL
        new_svlen_with_flag = 'SVLEN=' + '-' + str(2 * l_adj) #only for DEL
        info = info.replace(original_svlen_with_flag, new_svlen_with_flag)

    #Replace END value
    if info.split(';')[0].startswith('END='):
        #other info flags can contain "END=" ("CIEND=")
        original_end = abs(int(sv_line.split("\tEND=")[1].split(";")[0]))
        original_end_with_flag = "\tEND=" + str(original_end) + ";"
        converted_end_with_flag = "\tEND=" + str(converted_end) + ";"
        new_sv_line = "\t".join([region_id, str(converted_start), sv_id, ref, alt, qual, sv_filter, info]) 
        new_sv_line.replace(original_end_with_flag, converted_end_with_flag)
        return new_sv_line

    elif info.split(';')[-1].startswith('END='):
        original_end = abs(int(info.split("END=")[1]))
        original_end_with_flag = ";END=" + str(original_end)
        converted_end_with_flag = ";END=" + str(converted_end)

    else:
        original_end = abs(int(info.split(";END=")[1].split(";")[0]))
        original_end_with_flag = ";END=" + str(original_end) + ";"
        converted_end_with_flag = ";END=" + str(converted_end) + ";"

    new_info = info.replace(original_end_with_flag, converted_end_with_flag)
    new_sv_line = "\t".join([region_id, str(converted_start), sv_id, ref, alt, qual, sv_filter, new_info]) 
    return new_sv_line

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])