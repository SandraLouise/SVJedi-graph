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

    # parser.add_argument(
    #     "-l", 
    #     "--ladj", 
    #     metavar="<lengthOfFlankingRegions", 
    #     type=int,
    #     nargs=1, 
    #     default=[5000])
    
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
    # l_adj = int(args.ladj[0])


    # if args.outputDir:
    #     newRef_file_name = args.outputDir + '/reference_subsequences.fa'
    #     newVCF_file_name = args.outputDir + '/converted_variants.vcf'

    # else:
        ########### Tests
        # if args.number:
        #     n = int(args.number[0])
    
    d_breakpoints = {} 
    d_svs = {}
        
    if args.output:
        graph_file_name = args.output
        
    else:
        graph_file_name = inputVCF.split("/")[-1].replace(".vcf", "_graph.gfa")
        

    #1. Get reference genome sequence

    d_chrom = {}
    with open(reference_fasta) as sequenceFile:
        sequence = ""
        for line in sequenceFile:
            if line.startswith(">"):
                if sequence != "":
                    d_chrom[header] = sequence
                sequence = ""
                header = line.rstrip("\n")[1:].split()[0]

            else:
                sequence += line.rstrip("\n").upper()

        d_chrom[header] = sequence

    for chrom in d_chrom.keys():
        d_breakpoints[chrom] = set()
        d_svs[chrom] = []

    #2. Process data

    with open(inputVCF) as file:

        l_discarded = []
        dict_ins_seq = {}

        d_sv_alt = {}
        #d = {sv_id : [ref_seq, alt_seq]}

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

                # if not alt.startswith("<"):
                #     #ref and alt sequences provided in VCF line
                #     if ref[0] != alt[0]:
                #         start_on_chr -= 1
                #         dict_sv_alt[sv_id] = alt.upper()

                if sv_type == "DEL":
                    end_on_chr = int(get_info(info, "END"))
                    sv_id = format_DEL_INS_id(sv_type, pos, end_on_chr)
                
                elif sv_type == "INS":
                    #end_on_chr = int(get_info(info, "END"))
                    end_on_chr = start_on_chr
                    
                    sv_id = format_DEL_INS_id(sv_type, pos, end_on_chr)
                    #INS seq not in ALT field, get seq from INFO field
                    if alt.startswith("<"):
                        if any(["LEFT_SVINSSEQ=" in info, "RIGHT_SVINSSEQ=" in info]):
                            # left_insseq = get_info(info, "LEFT_SVINSSEQ")
                            # right_insseq = get_info(info, "RIGHT_SVINSSEQ")
                            l_discarded.append(line.rstrip())
                            continue

                        elif "SEQ=" in info:
                            dict_ins_seq[sv_id] = get_info(info, "SEQ")

                        else:
                            l_discarded.append(line.rstrip())
                            continue

                    #INS seq in ALT field but not in dict_sv_alt (no correction of pos required)
                    elif sv_id not in dict_ins_seq.keys():
                        dict_ins_seq[sv_id] = alt.upper()
                    
                elif sv_type == "INV":
                    end_on_chr = int(get_info(info, "END"))
                    sv_id = format_INV_id(pos, end_on_chr)
                
                elif sv_type == "BND":
                    alt = line.rstrip().split("\t")[4]
                    sv_id = format_BND_id(str(start_on_chr), alt)

                else: #SVTYPE = DUP par exemple
                    continue
                    print(sv_type)
    

                #Add breakpoints to list
                if sv_type in ['DEL', 'INS', 'INV']:

                    if end_on_chr >= len(d_chrom[chrom]) - 1 or start_on_chr >= len(d_chrom[chrom]) - 1:
                        l_discarded.append(line.rstrip())
                        continue

                    d_breakpoints[chrom].add(start_on_chr)
                    d_breakpoints[chrom].add(end_on_chr)

                    d_svs[chrom].append(sv_id)

                elif sv_type == 'BND':
                    left_coords, right_coords = parse_BND_id(chrom, sv_id)

                    if left_coords is None:
                        l_discarded.append(line.rstrip())
                        # continue

                    #Discard inter-chrom BNDs
                    elif left_coords[0] != right_coords[0]:
                        l_discarded.append(line.rstrip())
                        # continue

                    else:
                        if left_coords[2] == "-":
                            d_breakpoints[left_coords[0]].add(int(left_coords[1])-1)
                        else:
                            d_breakpoints[left_coords[0]].add(int(left_coords[1]))
                        if right_coords[2] == "+":
                            d_breakpoints[right_coords[0]].add(int(right_coords[1])-1)
                        else:
                            d_breakpoints[right_coords[0]].add(int(right_coords[1]))
                        
                        d_svs[chrom].append(sv_id)

                
                # if prev_chrom is None:
                    
                #     #Initialize breakpoints for new chrom
                #     region_breakpoints = [start_on_chr]
                #     if sv_type != 'BND':
                #         region_breakpoints.append(end_on_chr)

                #     #else: if 2nd bkpt BND sur meme chrom: ajout aux bkpt de la region

                #     region_id = '_'.join(['ref', chrom, sv_id])

                # elif prev_chrom != chrom:

                #     #Save breakpoints of prev chrom
                #     d_breakpoints[region_id] = region_breakpoints

                #     #Initialize breakpoints for new chrom
                #     region_breakpoints = [start_on_chr]
                #     if sv_type != 'BND':
                #         region_breakpoints.append(end_on_chr)

                #     #else: if 2nd bkpt BND sur meme chrom: ajout aux bkpt de la region

                #     region_id = '_'.join(['ref', chrom, sv_id])

                # else:
                #     region_breakpoints.append(start_on_chr)
                #     if sv_type != 'BND':
                #         region_breakpoints.append(end_on_chr)
                #     region_id = '_'.join([region_id, sv_id]) #ex: 'ref_1_DEL-88749-1400' + '_' + sv2_type + '-' + pos_sv_2 + '-' + len_sv_2

                # prev_chrom = chrom
        
        # # Last chrom
        # d_breakpoints[region_id] = region_breakpoints

        # Check if BND listed end in existing region (is there a region on endChr for which startRegion < endPos < endRegion ?)
        # if yes (expected case, translocations should be described with 2 BND lines): breakpoint should be in dict_breakpoints[region_id]
        # if not: add a region
    
    # List ignored SVs
    discarded = open('ignored_svs.vcf', 'w')
    discarded.write("##The following SVs were ignored during graph construction due to wrong format")
    for dis_sv in l_discarded:
        discarded.write("\n" + dis_sv)
    discarded.close()


    #3. Create graph
    graph_file = open(graph_file_name, "w")
    all_nodes = {}

    for chrom, sv_list in d_svs.items():
        if len(sv_list) > 0:
            graph_file.write("#{}\t{}\n".format(chrom, ";".join(sv_list)))

    for chrom, breakpoints in d_breakpoints.items():
        all_nodes[chrom] = []

        #Convert breakpoint set to breakpoint list + sort list
        breakpoints = list(breakpoints)
        breakpoints.sort()

        bkpt_to_remove = []
        for bkpt in breakpoints:
            if bkpt >= len(d_chrom[chrom]) -1:
                bkpt_to_remove.append(bkpt)
        for bkpt in bkpt_to_remove:
            breakpoints.remove(bkpt)

        bkpt_error = False

        # Write NODES, ref PATH and ref LINKS (excluding INS nodes and alternative links)
        nodes_len = []
        for i in range(-1, len(breakpoints)):

            # First node of region/graph component
            if i == -1:

                if len(breakpoints) == 0:
                    breakpoints.append(len(d_chrom[chrom]))

                #node_start = breakpoints[0]-l_adj
                node_start = 0
                node_end = breakpoints[0]-1

            # Last node of region/graph component
            elif i == len(breakpoints)-1:
                node_start = breakpoints[i]
                #node_end = breakpoints[i]+l_adj-1
                node_end = len(d_chrom[chrom]) - 1

            # INS node added later
            # elif breakpoints[i] == breakpoints[i+1]:
            #     continue

            else:
                node_start = breakpoints[i]
                node_end = breakpoints[i+1]-1
            
            node_id = format_node_id(chrom, node_start, node_end)
            all_nodes[chrom].append(node_id)
            node_seq = get_ref_seq(d_chrom, chrom, node_start, node_end)
            if node_end - node_start + 1 != len(node_seq):
                print("Error in node:", chrom, str(node_start), str(node_end), str(len(node_seq)))
                bkpt_error = True

            graph_file.write(format_gfa_node(node_id, node_seq))
            nodes_len.append(str(len(node_seq)))

        path_nodes = "+,".join(all_nodes[chrom]) + "+"
        path_len = "M,".join(nodes_len) + "M"
        graph_file.write(format_gfa_path(chrom, path_nodes, path_len))

        if bkpt_error:
            print(chrom)
            print(breakpoints)

        # print(region, all_nodes[chrom])
        # print(path_nodes)

        for i in range(len(all_nodes[chrom])-1):
            # print(all_nodes[chrom][i], all_nodes[chrom][i+1])
            graph_file.write(format_gfa_link_pp(all_nodes[chrom][i], all_nodes[chrom][i+1], 0))
        
        # all_nodes.extend(all_nodes[chrom])

    #Write alt NODES and alt LINKS
    for chrom in d_svs.keys():
        svs = d_svs[chrom]
        for sv_id in svs:
            sv_type = sv_id.split("-")[0]

            if sv_type != "BND":
                pos, end = sv_id.split("-")[1:]
            else:
                alt = sv_id.split("-")[1]

            pos = int(pos)

            if sv_type == "DEL":
                # end = pos + int(len_end) # /!\ end pos of nodes is defined from END tag in info field of vcf !!
                
                # for node in all_nodes[chrom]:
                #     coords = node.split(":")[1]
                #     if any([coords.startswith(str(pos)), coords.startswith(str(pos+1))]):
                #         end = int(coords.split("-")[1])
                #         # => PB if another breakpoint at pos+1 !!

                if sv_id in d_sv_alt.keys():
                    #alternative node
                    alt_node = format_node_id(chrom, pos-1, "a")
                    graph_file.write(format_gfa_node(alt_node, d_sv_alt[sv_id]))
                    #alternative links
                    for node in all_nodes[chrom]:
                        coords = node.split(":")[1]
                        if any([coords.endswith(str(pos-1)), coords.endswith(str(pos))]):
                            left_node = node
                        elif coords.startswith(str(end)):
                            right_node = node
                    graph_file.write(format_gfa_link_pp(left_node, alt_node, 1))
                    graph_file.write(format_gfa_link_pp(alt_node, right_node, 1))

                else:
                    end = int(end)
                    for node in all_nodes[chrom]:
                        coords = node.split(":")[1]
                        # if coords.endswith(str(pos)) or coords.endswith(str(pos-1)):
                        if coords.endswith(str(pos)):
                            left_node = node
                        elif coords.startswith(str(end+1)):
                            right_node = node
                    graph_file.write(format_gfa_link_pp(left_node, right_node, 1))
            
            elif sv_type == "INS":

                left_node, right_node = None, None

                if sv_id in d_sv_alt.keys():
                    # corr_pos = pos - 1
                    ins_node = format_altnode_id(chrom, pos+1)
                    graph_file.write(format_gfa_node(ins_node, d_sv_alt[sv_id]))

                    for node in all_nodes[chrom]:
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
                    graph_file.write(format_gfa_node(ins_node, dict_ins_seq[sv_id]))

                    for node in all_nodes[chrom]:
                        coords = node.split(":")[1]
                        # if any([coords.endswith(str(pos)), coords.endswith(str(pos-1))]):
                        if coords.endswith(str(pos)):
                            left_node = node
                        # elif any([coords.startswith(str(pos+1)), coords.startswith(str(pos))]):
                        elif coords.startswith(str(pos+1)):
                            right_node = node

                    # if any([left_node is None, right_node is None]):
                        # print(sv, all_nodes[chrom])
                    graph_file.write(format_gfa_link_pp(left_node, ins_node, 1))
                    graph_file.write(format_gfa_link_pp(ins_node, right_node, 1))

            elif sv_type == "INV":
                end = int(end)

                left_node, left_invnode, right_invnode, right_node = [None]*4
                for node in all_nodes[chrom]:
                    start_pos, end_pos = node.split(":")[1].split("-")
                    start_pos, end_pos = int(start_pos), int(end_pos)
                    if end_pos == pos:
                        left_node = node
                    elif start_pos == end + 1:
                        right_node = node
                    else:
                        if start_pos == pos + 1:
                            left_invnode = node
                        if end_pos == end:
                            right_invnode = node
                
                if any([left_node is None, right_node is None, left_invnode is None, right_invnode is None]):
                    # print(sv, "\tleft:", left_node, "\tinv:", left_invnode, right_invnode, "\tright:", right_node)
                    continue

                else:
                    # print(sv, "\tleft:", left_node, "\tinv:", left_invnode, right_invnode, "\tright:", right_node)
                    graph_file.write(format_gfa_link_pm(left_node, right_invnode))
                    graph_file.write(format_gfa_link_mp(left_invnode, right_node))
            
            elif sv_type == "BND":
                left_coords, right_coords = parse_BND_id(chrom, sv_id)

                if left_coords[2] == "-":
                    left_node = find_node_by_start(left_coords[0], int(left_coords[1]), all_nodes[left_coords[0]])
                else:
                    left_node = find_node_by_end(left_coords[0], int(left_coords[1]), all_nodes[left_coords[0]])

                if right_coords[2] == "+":
                    right_node = find_node_by_start(right_coords[0], int(right_coords[1]), all_nodes[right_coords[0]])
                else:
                    right_node = find_node_by_end(right_coords[0], int(right_coords[1]), all_nodes[right_coords[0]])
                
                if left_coords[2] == "-":
                    graph_file.write(format_gfa_link_mp(left_node, right_node))
                elif right_coords[2] == "-":
                    graph_file.write(format_gfa_link_pm(left_node, right_node))
                else:
                    graph_file.write(format_gfa_link_pp(left_node, right_node, 1))

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

def find_node_by_start(Tchr, pos, node_list):
    for node in node_list:
        chrom, coords = node.split(":")
        start_pos = int(coords.split("-")[0])
        if chrom == Tchr and start_pos == pos:
            return node
    #Test
    sys.exit("Node not found, incorrect start position: " + str(pos))

def find_node_by_end(Tchr, pos, node_list):
    for node in node_list:
        chrom, coords = node.split(":")
        end_pos = int(coords.split("-")[1])
        if chrom == Tchr and end_pos == pos:
            return node
    #Test
    sys.exit("Node not found, incorrect end position: " + str(pos))

def format_DEL_INS_id(sv_type, pos, end):
    return '-'.join([sv_type, str(pos), str(end)])
def format_INV_id(pos, end):
    return '-'.join(["INV", str(pos), str(end)])
def format_BND_id(pos, alt):

    if "[" in alt:
        parsed_alt = list(filter(bool, alt.split("[")))
        if ":" in parsed_alt[1]:
            s = parsed_alt[0]
        else:
            s = parsed_alt[1]

    elif "]" in alt:
        parsed_alt = list(filter(bool, alt.split("]")))
        if ":" in parsed_alt[1]:
            s = parsed_alt[0]
        else:
            s = parsed_alt[1]

    else:
        return "BND-format"

    if len(s) > 1:
        return "BND-format"

    alt = alt.replace(s, pos)
    return '-'.join(["BND", alt])

def parse_BND_id(chrom, bnd_id):

    right_chrom = None
    alt = bnd_id.split("BND-")[1]

    if "[" in alt:
        alt = list(filter(bool, alt.split("[")))
        if ":" in alt[1]:
            # t[p[ : piece extending to the right of p is joined after t
            right_chrom = alt[1].split(":")[0]
            right_pos = alt[1].split(":")[1]
            right_strand = "+"
            left_chrom = chrom
            left_pos = alt[0]
            left_strand = "+"

        elif ":" in alt[0]:
            # [p[t : reverse comp piece extending right of p is joined before t
            left_chrom = alt[0].split(":")[0]
            left_pos = alt[0].split(":")[1]
            left_strand = "-"
            right_chrom = chrom
            right_pos = alt[1]
            right_strand = "+"

    elif "]" in alt:
        alt = list(filter(bool, alt.split("]")))
        if ":" in alt[1]:
            # t]p] : reverse comp piece extending left of p is joined after t
            right_chrom = alt[1].split(":")[0]
            right_pos = alt[1].split(":")[1]
            right_strand = "-"
            left_chrom = chrom
            left_pos = alt[0]
            left_strand = "+"

        elif ":" in alt[0]:
            # ]p]t piece extending to the left of p is joined before t
            left_chrom = alt[0].split(":")[0]
            left_pos = alt[0].split(":")[1]
            left_strand = "+"
            right_chrom = chrom
            right_pos = alt[1]
            right_strand = "+"

    if right_chrom is None:
        return None, None
    else:
        return (left_chrom, left_pos, left_strand), (right_chrom, right_pos, right_strand)


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