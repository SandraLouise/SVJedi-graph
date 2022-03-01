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
        default=[5000])
    
    parser.add_argument(
        "-o", 
        "--output", 
        metavar="<outputFile", 
        type=str)
    
    # parser.add_argument(
    #     "-o", 
    #     "--outputDir", 
    #     metavar="<outputDirectory>", 
    #     type=str,
    #     required=False)

    args = parser.parse_args()

    inputVCF = args.vcf[0]
    reference_fasta = args.ref[0]
    l_adj = int(args.ladj[0])


    # if args.outputDir:
    #     newRef_file_name = args.outputDir + '/reference_subsequences.fa'
    #     newVCF_file_name = args.outputDir + '/converted_variants.vcf'

    # else:
        ########### Tests
        # if args.number:
        #     n = int(args.number[0])
        
    if args.output:
        graph_file_name = args.output
        
    else:
        graph_file_name = inputVCF.split("/")[-1].replace(".vcf", "_graph.gfa")


    #1. Get reference genome sequence

    dict_of_chrom = {}
    with open(reference_fasta) as sequenceFile:
        sequence = ""
        for line in sequenceFile:
            if line.startswith(">"):
                if sequence != "":
                    dict_of_chrom[header] = sequence
                sequence = ""
                header = line.rstrip("\n")[1:].split(" ")[0]

            else:
                sequence += line.rstrip("\n").upper()

        dict_of_chrom[header] = sequence


    #2. Process data

    with open(inputVCF) as file:

        dict_breakpoints = {} 
        #d = {id_region : [svStart_onChr, svEnd_onChr, ...]}

        dict_sv_alt = {}
        #d = {sv_id : [ref_seq, alt_seq]}

        dict_ins_seq = {}
        bnd_list = []
        #format: startChr:startPos;alt

        prev_chrom = None

        for line in file:
            if line.startswith('#'):
                continue

            else:
                #Test
                # if args.number:
                #     if count == n:
                #         break

                #Get SV info
                chrom, pos, sv_id, ref, alt, __, __, info, *__ = line.rstrip().split("\t")
                sv_type = get_info(info, "SVTYPE")
                start_on_chr = int(pos) 
                vcf_id = sv_id

                if not alt.startswith("<"):
                    #ref and alt sequences provided in VCF line
                    if ref[0] != alt[0]:
                        start_on_chr -= 1
                        dict_sv_alt[sv_id] = alt.upper()

                if sv_type == "DEL":
                    end_on_chr = int(get_info(info, "END"))
                    sv_id = format_DEL_INS_id(sv_type, pos, end_on_chr)
                
                elif sv_type == "INS":
                    #end_on_chr = int(get_info(info, "END"))
                    end_on_chr = start_on_chr
                    
                    sv_id = format_DEL_INS_id(sv_type, pos, end_on_chr)
                    #INS seq not in ALT field, get seq from INFO field
                    if alt.startswith("<"):
                        dict_ins_seq[sv_id] = get_info(info, "SEQ")
                    #INS seq in ALT field but not in dict_sv_alt (no correction of pos required)
                    elif sv_id not in dict_sv_alt.keys():
                        dict_ins_seq[sv_id] = alt.upper()
                    
                elif sv_type == "INV":
                    end_on_chr = int(info.split(';END=')[1].split(';')[0])
                    sv_id = format_INV_id(pos, end_on_chr)
                
                elif sv_type == "BND":
                    continue

                    # alt = line.rstrip().split("\t")[4]
                    # bnd_list.append(";".join([":".join([chrom, str(pos)]), alt])) #chrom:pos;alt

                    # if "[" in alt:
                    #     end_chrom = alt.split(":")[0].split("[")[1]
                    #     end_pos = alt.split(":")[1].split("[")[0]
                    # elif "]" in alt:
                    #     end_chrom = alt.split(":")[0].split("]")[1]
                    #     end_pos = alt.split(":")[1].split("]")[0]
                    
                    # sv_id = format_BND_id(pos, end_chrom, end_pos)

                else: #SVTYPE = DUP par exemple
                    continue
                
                if prev_chrom is None:
                    
                    #Initialize breakpoints for new chrom
                    region_breakpoints = [start_on_chr]
                    if sv_type != 'BND':
                        region_breakpoints.append(end_on_chr)

                    #else: if 2nd bkpt BND sur meme chrom: ajout aux bkpt de la region

                    region_id = '_'.join(['ref', chrom, sv_id])

                elif prev_chrom != chrom:

                    #Save breakpoints of prev chrom
                    dict_breakpoints[region_id] = region_breakpoints

                    #Initialize breakpoints for new chrom
                    region_breakpoints = [start_on_chr]
                    if sv_type != 'BND':
                        region_breakpoints.append(end_on_chr)

                    #else: if 2nd bkpt BND sur meme chrom: ajout aux bkpt de la region

                    region_id = '_'.join(['ref', chrom, sv_id])

                else:
                    region_breakpoints.append(start_on_chr)
                    if sv_type != 'BND':
                        region_breakpoints.append(end_on_chr)
                    region_id = '_'.join([region_id, sv_id]) #ex: 'ref_1_DEL-88749-1400' + '_' + sv2_type + '-' + pos_sv_2 + '-' + len_sv_2

                prev_chrom = chrom
        
        # Last chrom
        dict_breakpoints[region_id] = region_breakpoints

        # Check if BND listed end in existing region (is there a region on endChr for which startRegion < endPos < endRegion ?)
        # if yes (expected case, translocations should be described with 2 BND lines): breakpoint should be in dict_breakpoints[region_id]
        # if not: add a region


    #3. Create graph
    graph_file = open(graph_file_name, "w")
    all_nodes = []

    for region, breakpoints in dict_breakpoints.items():
        chrom = region.split("_")[1]
        svs = region.split("_")[2:]

        breakpoints.sort()

        # Write NODES, ref PATH and ref LINKS (excluding INS nodes and alternative links)
        ref_nodes = []
        nodes_len = []
        for i in range(-1, len(breakpoints)):
            if i == -1:
                node_start = breakpoints[0]-l_adj
                node_end = breakpoints[0]-1
            elif i == len(breakpoints)-1:
                node_start = breakpoints[i]
                node_end = breakpoints[i]+l_adj-1
            elif breakpoints[i] == breakpoints[i+1]: # skip INS node
                continue
            else:
                node_start = breakpoints[i]
                node_end = breakpoints[i+1]-1
            
            node_id = format_node_id(chrom, node_start, node_end)
            ref_nodes.append(node_id)
            node_seq = get_ref_seq(dict_of_chrom, chrom, node_start, node_end)
            if node_end - node_start + 1 != len(node_seq):
                print("Error in node:", chrom, str(node_start), str(node_end), str(len(node_seq)))

            graph_file.write(format_gfa_node(node_id, node_seq))
            nodes_len.append(str(len(node_seq)))

        path_nodes = "+,".join(ref_nodes) + "+"
        path_len = "M,".join(nodes_len) + "M"
        graph_file.write(format_gfa_path(region, path_nodes, path_len))

        # print(region, ref_nodes)
        # print(path_nodes)

        for i in range(len(ref_nodes)-1):
            # print(ref_nodes[i], ref_nodes[i+1])
            graph_file.write(format_gfa_link_pp(ref_nodes[i], ref_nodes[i+1], 0))
        
        all_nodes.extend(ref_nodes)

        #Write alt NODES and alt LINKS
        for sv in svs:
            sv_type, pos, end = sv.split("-")
            pos = int(pos)

            if sv_type == "DEL":
                # end = pos + int(len_end) # /!\ end pos of nodes is defined from EnD tag in info field of vcf !!
                for node in ref_nodes:
                    coords = node.split(":")[1]
                    if any([coords.startswith(str(pos)), coords.startswith(str(pos+1))]):
                        end = int(coords.split("-")[1])

                if sv in dict_sv_alt.keys():
                    #alternative node
                    alt_node = format_node_id(chrom, pos-1, "a")
                    graph_file.write(format_gfa_node(alt_node, dict_sv_alt[sv]))
                    #alternative links
                    for node in ref_nodes:
                        coords = node.split(":")[1]
                        if any([coords.endswith(str(pos-1)), coords.endswith(str(pos))]):
                            left_node = node
                        elif coords.startswith(str(end)):
                            right_node = node
                    graph_file.write(format_gfa_link_pp(left_node, alt_node, 1))
                    graph_file.write(format_gfa_link_pp(alt_node, right_node, 1))

                else:
                    for node in ref_nodes:
                        coords = node.split(":")[1]
                        if coords.endswith(str(pos)) or coords.endswith(str(pos-1)):
                            left_node = node
                        elif coords.startswith(str(end+1)):
                            right_node = node
                    graph_file.write(format_gfa_link_pp(left_node, right_node, 1))
            
            elif sv_type == "INS":

                left_node, right_node = None, None

                if sv in dict_sv_alt.keys():
                    # corr_pos = pos - 1
                    ins_node = format_altnode_id(chrom, pos+1)
                    graph_file.write(format_gfa_node(ins_node, dict_sv_alt[sv]))

                    for node in ref_nodes:
                        coords = node.split(":")[1]
                        if coords.endswith(str(pos)):
                            left_node = node
                        elif coords.startswith(str(pos+1)):
                            right_node = node
                    graph_file.write(format_gfa_link_pp(left_node, ins_node, 1))
                    graph_file.write(format_gfa_link_pp(ins_node, right_node, 1))

                else:
                    # Coordinates of ins sequence end with "a" (for alternative)
                    ins_node = format_altnode_id(chrom, pos+1)
                    graph_file.write(format_gfa_node(ins_node, dict_ins_seq[sv]))

                    for node in ref_nodes:
                        coords = node.split(":")[1]
                        if any([coords.endswith(str(pos)), coords.endswith(str(pos-1))]):
                            left_node = node
                        elif any([coords.startswith(str(pos+1)), coords.startswith(str(pos))]):
                            right_node = node

                    # if any([left_node is None, right_node is None]):
                        # print(sv, ref_nodes)
                    graph_file.write(format_gfa_link_pp(left_node, ins_node, 1))
                    graph_file.write(format_gfa_link_pp(ins_node, right_node, 1))

            elif sv_type == "INV":
                end = int(end)
                for node in ref_nodes:
                    coords = node.split(":")[1]
                    if coords.endswith(str(pos-1)):
                        left_node = node
                    elif coords.startswith(str(end)):
                        right_node = node
                    elif coords.startswith(str(pos)):
                        inv_node = node
                graph_file.write(format_gfa_link_pm(left_node, inv_node))
                graph_file.write(format_gfa_link_mp(inv_node, right_node))
        
    # Write alt links for BNDs
    for bnd in bnd_list:
        chrom, pos = bnd.split(";")[0].split(":")
        pos = int(pos)
        alt = bnd.split(";")[1]

        if "]" in alt:
            end_chr = alt.split(":")[0].split("]")[1]
            end_pos = int(alt.split(":")[1].split("]")[0])
        elif "[" in alt:
            end_chr = alt.split(":")[0].split("[")[1]
            end_pos = int(alt.split(":")[1].split("[")[0])

        if alt.endswith('['):
            #t[p[
            #piece extending to the right of p is joined after t
            left_node = find_node_by_end(chrom, pos, all_nodes) #node ending at chr1:pos
            right_node = find_node_by_start(end_chr, end_pos, all_nodes) #node starting at end_chr:end_pos
            graph_file.write(format_gfa_link_pp(left_node, right_node))
                    
        elif alt.endswith(']'):
            #t]p]   
            #reverse comp piece extending left of p is joined after t   
            left_node = find_node_by_end(chrom, pos-1, all_nodes) #node ending at chr1:pos
            right_node = find_node_by_end(end_chr, end_pos-1, all_nodes) #rev of node ending at end_chr:end_pos -1
            graph_file.write(format_gfa_link_pm(left_node, right_node))
                        
        elif alt.startswith(']'):
            #]p]t
            #piece extending to the left of p is joined before t 
            left_node = find_node_by_end(end_chr, end_pos, all_nodes) #node ending at end_chr:end_pos
            right_node = find_node_by_start(chrom, pos, all_nodes) #node starting at chr1:pos
            graph_file.write(format_gfa_link_pp(left_node, right_node))
                    
        elif alt.startswith('['):
            #[p[t[
            #reverse comp piece extending right of p is joined before t 
            left_node = find_node_by_start(end_chr, end_pos, all_nodes) #rev of node starting at end_chr:end_pos
            right_node = find_node_by_start(chrom, pos, all_nodes) #node starting at chr1:pos
            graph_file.write(format_gfa_link_mp(left_node, right_node))


    # Add BND links
    # for bnd in bnd_list:   
    #   add 1 link
    #   (find node_id: (1) find region with first_sv_pos < bnd_pos < last_sv_pos (2) from region ref path, find node_id starting/ending with bnd_pos)    

    graph_file.close()


def format_gfa_link_pp(node1, node2, rank):
    # print("\t".join(["L", node1, "+", node2, "+", "0M", "SR:i:"+str(rank)]) + "\n")
    return "\t".join(["L", node1, "+", node2, "+", "0M"]) + "\n"
def format_gfa_link_pm(node1, node2):
    return "\t".join(["L", node1, "+", node2, "-", "0M"]) + "\n"
def format_gfa_link_mp(node1, node2):
    return "\t".join(["L", node1, "-", node2, "+", "0M"]) + "\n"
def format_gfa_link_mm(node1, node2):
    return "\t".join(["L", node1, "-", node2, "-", "0M"]) + "\n"

def format_gfa_node(node_id, seq):
    '''Convert node info in GFA format. For node_id format pattern, see function format_node_id().'''
    return "\t".join(["S", node_id, seq]) + "\n"
def format_gfa_path(name, nodes, lengths):
    return "\t".join(["P", name, nodes, lengths]) + "\n"

def format_node_id(chrom, coord1, coord2):
    # Add +1 to start and end positions to use 1-indexed positions (compatible with 1-indexed positions of VCF file format)
    return ":".join([ str(chrom), "-".join([ str(coord1+1), str(coord2+1) ]) ])
def format_altnode_id(chrom, coord1):
    return ":".join([ str(chrom), str(coord1)+"a" ])

def find_node_by_start(chr, pos, node_list):
    for node in node_list:
        chrom, coords = node.split(":")
        if chrom == chr and coords.startswith(str(pos)):
            return node
    #Test
    sys.exit("Node not found, incorrect start position given: " + str(pos))

def find_node_by_end(Tchr, pos, node_list):
    for node in node_list:
        chrom, coords = node.split(":")
        if chrom == Tchr and coords.endswith(str(pos)):
            return node
    #Test
    sys.exit("Node not found, incorrect end position given: " + str(pos))

def format_DEL_INS_id(sv_type, pos, end):
    return '-'.join([sv_type, str(pos), str(end)])
def format_INV_id(pos, end):
    return '-'.join(["INV", str(pos), str(end)])
def format_BND_id(pos, end_chr, end_pos):
    return '-'.join(["BND", str(pos), ":".join([end_chr, end_pos])])

def get_info(info, label):
    #Info label first in info field
    if info.split(";")[0].startswith(label+"="):
        return info.split(label+"=")[1].split(";")[0]

    #Info label last in info field
    elif info.split(";")[-1].startswith(label+"="):
        return info.split(label.join([";", "="]))[1]

    else:
        return info.split(label.join([";", "="]))[1].split(";")[0]

def get_sv_len(sv_info):
    if sv_info.split(';')[-1].startswith('SVLEN='):
        return abs(int(sv_info.split("SVLEN=")[1]))
    else:
        return abs(int(sv_info.split("SVLEN=")[1].split(";")[0]))

# def get_sv_type(sv_info):
#     if 'SVTYPE' in sv_info:
#         if sv_info.split(';')[-1].startswith('SVTYPE='):
#             return sv_info.split('SVTYPE=')[1]
#         else:
#             return sv_info.split('SVTYPE=')[1].split(';')[0]

def get_ref_seq(dict_of_chrom, chrom, start_included, end_included):
    return dict_of_chrom[chrom][start_included : end_included + 1]


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])
