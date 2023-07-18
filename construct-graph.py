#!/usr/bin/env python3

"""*******************************************************************************
    Name: SVJedi-graph
    Description: SVjedi-graph aims to genotype structural variant with long reads data using a variation graph.
    Author: Sandra Romain
    Contact: sandra.romain@inria.fr, Inria/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
    
    Copyright (C) 2022 Inria
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""

import sys
import argparse
import json

def parse_arguments(arguments):

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
        "-o", 
        "--output", 
        metavar="<outputFile", 
        type=str)
    
    args = parser.parse_args(arguments)

    inVCF = args.vcf[0]
    inFA = args.ref[0]

    if args.output:
        outGFA = args.output
        
    else:
        outGFA = inVCF.split("/")[-1].replace(".vcf", "_graph.gfa")

    return inVCF, inFA, outGFA

def construct_gfa(inVCF, inFA, outGFA):

    d_chr_bkpt = {} # (key) chromosome -> (value) list of breakpoint positions (int)
    d_bkpt_sv = {} # (key) chromosome -> (value) { (key) breakpoint position -> (value) list of sv_id }
    d_link_sv = {} # (key) link (node1, orient1, node2, orient2) -> (value) list of sv_id
    d_svs = {}
    d_sv_ID = {} # (key) sv_id -> (value) vcf_id
    
    #================================================================================
    # 1. Get reference genome sequence
    #================================================================================

    d_chrom = {}
    with open(inFA) as sequenceFile:
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
        d_chr_bkpt[chrom] = set()
        d_svs[chrom] = []

    #================================================================================
    # 2. Process data
    #================================================================================

    with open(inVCF) as file:

        l_discarded = []
        dict_ins_seq = {}
        d_ins_multiplicity = {}

        d_sv_alt = {}
        #d = {sv_id : [ref_seq, alt_seq]}

        for line in file:
            if line.startswith('#'):
                continue

            else:

                #--------------------------------------------------------------------
                #Get SV info
                #--------------------------------------------------------------------
                chrom, pos, sv_id, ref, alt, __, __, info, *__ = line.rstrip().split("\t")
                sv_type = get_info(info, "SVTYPE")
                start_on_chr = int(pos) 
                vcf_id = sv_id

                if sv_type == "DEL":
                    end_on_chr = int(get_info(info, "END"))
                    sv_id = format_DEL_id(pos, end_on_chr)
                
                elif sv_type == "INS":
                    end_on_chr = start_on_chr
                    
                    # Handle multiple INS at same position
                    if pos not in d_ins_multiplicity.keys():
                        d_ins_multiplicity[pos] = 0
                    d_ins_multiplicity[pos] += 1
                    ins_count = d_ins_multiplicity[pos]
                    
                    sv_id = format_INS_id(pos, ins_count)

                    #Wrong format: REF field len > 1
                    if len(ref) > 1 :
                        l_discarded.append(line.rstrip())
                        continue

                    #INS seq not in ALT field, get seq from INFO field
                    elif alt.startswith("<"):
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

                d_sv_ID[sv_id] = vcf_id

                ## DEBUG
                # if sv_type == "BND":
                #     for coord in parse_BND_id(chrom, sv_id):
                #         if coord[1] == 1:
                #             print(line)
                # else:
                #     if 1 in [start_on_chr, end_on_chr]:
                #         print(line)
                ##
    
                #--------------------------------------------------------------------
                #Add breakpoints to list
                #--------------------------------------------------------------------
                if sv_type in ['DEL', 'INS', 'INV']:

                    if end_on_chr >= len(d_chrom[chrom]) - 1 or start_on_chr >= len(d_chrom[chrom]) - 1:
                        l_discarded.append(line.rstrip())
                        continue

                    # Add chrom to d_bkpt_sv keys
                    if chrom not in d_bkpt_sv.keys():
                        d_bkpt_sv[chrom] = {}

                    for bkpt_pos in set([start_on_chr, end_on_chr]):
                        
                        # Do not add bkpt pos at extremities of chrom
                        if 1 < bkpt_pos < len(d_chrom[chrom]):

                            # Add breakpoint position to chromosome breakpoint list
                            d_chr_bkpt[chrom].add(bkpt_pos)

                            # Add breakpoint position to d_bkpt_sv[chrom] keys
                            if bkpt_pos not in d_bkpt_sv[chrom].keys():
                                d_bkpt_sv[chrom][bkpt_pos] = []

                            # Associate breakpoint position to SV
                            d_bkpt_sv[chrom][bkpt_pos].append(sv_id)

                    d_svs[chrom].append(sv_id)

                elif sv_type == 'BND':
                    left_coords, right_coords = parse_BND_id(chrom, sv_id)

                    # Unsupported BND format
                    if left_coords is None:
                        l_discarded.append(line.rstrip())
                        # continue

                    else:
                        #--------------------------------------------------------------------
                        # Correct bkpt pos
                        #--------------------------------------------------------------------
                        # Case 1: left and right are forward strands (t[p[ or ]p]t)
                        if left_coords[2] == "+" and right_coords[2] == "+":

                            # The bkpt POS on right chrom is p -1
                            right_coords[1] = right_coords[1] - 1
                        
                        # Case 2: right is reverse strand (t]p])
                            # p is the bkpt POS on right chrom
                            # no need to modify right_coords[1]
                        
                        # Case 3: left is reverse strand ([p[t)
                        elif left_coords[2] == "-":

                            # The bkpt POS on left chrom is p - 1
                            left_coords[1] = left_coords[1] - 1
                            # The bkpt POS on right chrom is t - 1
                            right_coords[1] = right_coords[1] - 1

                        # Get the 2 breakpoints
                        left_bkpt = (left_coords[0], left_coords[1])
                        right_bkpt = (right_coords[0], right_coords[1])

                        # Add breakpoints to dictionnaries
                        for (bkpt_chrom, bkpt_pos) in [left_bkpt, right_bkpt]:

                            # Do not add bkpt pos at extremities of chrom
                            if 1 < bkpt_pos < len(d_chrom[bkpt_chrom]):
                            
                                # Add breakpoint position to chromosome breakpoint list
                                if bkpt_chrom not in d_bkpt_sv.keys():
                                    d_bkpt_sv[bkpt_chrom] = {}

                                d_chr_bkpt[bkpt_chrom].add(bkpt_pos)

                                # Add breakpoint position to d_bkpt_sv[chrom] keys
                                if bkpt_pos not in d_bkpt_sv[bkpt_chrom].keys():
                                    d_bkpt_sv[bkpt_chrom][bkpt_pos] = []
                                
                                # Associate breakpoint position to SV
                                d_bkpt_sv[bkpt_chrom][bkpt_pos].append(sv_id)

                        d_svs[chrom].append(sv_id)
    
    # List ignored SVs
    discarded = open('ignored_svs.vcf', 'w')
    discarded.write("##The following SVs were ignored during graph construction due to wrong format")
    for dis_sv in l_discarded:
        discarded.write("\n" + dis_sv)
    discarded.close()

    #================================================================================
    # 3. Create graph
    #================================================================================

    graph_file = open(outGFA, "w")
    all_nodes = {}

    for chrom, sv_list in d_svs.items():
        if len(sv_list) > 0:
            graph_file.write("#{}\t{}\n".format(chrom, ";".join(sv_list)))

    for chrom, breakpoints in d_chr_bkpt.items():
        all_nodes[chrom] = []

        #----------------------------------------------------------------------------
        # Convert breakpoint set to breakpoint list + sort list
        #----------------------------------------------------------------------------
        breakpoints = list(breakpoints)
        breakpoints.sort()

        bkpt_to_remove = []
        for bkpt in breakpoints:
            if bkpt >= len(d_chrom[chrom]) -1:
                bkpt_to_remove.append(bkpt)
        for bkpt in bkpt_to_remove:
            breakpoints.remove(bkpt)

        bkpt_error = False

        #----------------------------------------------------------------------------
        # Write NODES, ref PATH and ref LINKS (excluding INS nodes and alternative links)
        #----------------------------------------------------------------------------
        nodes_len = []
        for i in range(-1, len(breakpoints)):

            # First node of region/graph component
            if i == -1:

                if len(breakpoints) == 0:
                    breakpoints.append(len(d_chrom[chrom]))

                node_start = 0
                node_end = breakpoints[0]-1

            # Last node of region/graph component
            elif i == len(breakpoints)-1:
                node_start = breakpoints[i]
                node_end = len(d_chrom[chrom]) - 1

            # INS node added later
            # elif breakpoints[i] == breakpoints[i+1]:
            #     continue

            else:
                node_start = breakpoints[i]
                node_end = breakpoints[i+1]-1

            # Create reference node
            node_id = format_node_id(chrom, node_start, node_end)
            all_nodes[chrom].append(node_id)
            node_seq = get_ref_seq(d_chrom, chrom, node_start, node_end)

            # Check for inaccuracies between node len and node positions
            if node_end - node_start + 1 != len(node_seq):
                print(f"Warning: seq. length of node {node_id} ({str(len(node_seq))}) doesn't match node coordinates ({str(node_start)} to {str(node_end)})")
                bkpt_error = True

            # Add node to graph (gfa)
            graph_file.write(format_gfa_node(node_id, node_seq))
            nodes_len.append(str(len(node_seq)))

            # Add link to d_link_sv
            if i != -1:

                # Create reference link (between previous node and current node) and add to graph (gfa)
                graph_file.write(format_gfa_link_pp(prev_node_id, node_id, 0))

                # Associate link to sv_id (reference allele)
                link_id = (prev_node_id, "+", node_id, "+")

                link_key = get_link_key(link_id)

                d_link_sv[link_key] = []

                for sv_id in d_bkpt_sv[chrom][breakpoints[i]]:
                    d_link_sv[link_key].append((":".join([chrom, sv_id]), 0))
            
            # Save node_id as prev_node_id
            prev_node_id = node_id

        path_nodes = "+,".join(all_nodes[chrom]) + "+"
        path_len = "M,".join(nodes_len) + "M"
        graph_file.write(format_gfa_path(chrom, path_nodes, path_len))

        if bkpt_error:
            print(chrom)
            print(breakpoints)
    
    #----------------------------------------------------------------------------
    # Write alt NODES and alt LINKS
    #----------------------------------------------------------------------------
    for chrom, svs in d_svs.items():

        for sv_id in svs:
            sv_type = sv_id.split("-")[0]

            if sv_type == "INS":
                pos, ins_count = sv_id.split("-")[1:]
                pos = int(pos)

            elif sv_type == "BND":
                alt = sv_id.split("-")[1]
                
            else:
                pos, end = sv_id.split("-")[1:]
                pos, end = int(pos), int(end)


            if sv_type == "DEL":

                #--------------------------------------------------------------------------
                # Part to remove
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
                #--------------------------------------------------------------------------

                else:

                    for node in all_nodes[chrom]:
                        coords = node.split(":")[1]
                        # if coords.endswith(str(pos)) or coords.endswith(str(pos-1)):
                        if coords.endswith(str(pos)):
                            left_node = node
                        elif coords.startswith(str(end+1)):
                            right_node = node
                    
                    # Create alternative link and add to graph (gfa)
                    graph_file.write(format_gfa_link_pp(left_node, right_node, 1))

                    # Associate alternative link to sv_id (DEL)
                    link_id = (left_node, "+", right_node, "+")

                    link_key = get_link_key(link_id)

                    if link_key not in d_link_sv.keys():
                        d_link_sv[link_key] = []
                    
                    d_link_sv[link_key].append((":".join([chrom, sv_id]), 1))
            
            elif sv_type == "INS":

                left_node, right_node = None, None

                #--------------------------------------------------------------------------
                # Part to remove
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
                #--------------------------------------------------------------------------

                else:
                    # Coordinates of ins sequence end with "a" (for alternative)
                    ins_node = format_altnode_id(chrom, pos+1, ins_count)
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

                    # Create alternative links and add to graph (gfa)
                    graph_file.write(format_gfa_link_pp(left_node, ins_node, 1))
                    graph_file.write(format_gfa_link_pp(ins_node, right_node, 1))

                    # Associate alternative links to sv_id (INS)
                    link_id_1 = (left_node, "+", ins_node, "+")
                    link_id_2 = (ins_node, "+", right_node, "+")

                    for link_id in [link_id_1, link_id_2]:

                        link_key = get_link_key(link_id)

                        if link_key not in d_link_sv.keys():
                            d_link_sv[link_key] = []
                        
                        d_link_sv[link_key].append((":".join([chrom, sv_id]), 1))

            elif sv_type == "INV":

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

                    # Create alternative links and add to graph (gfa)
                    graph_file.write(format_gfa_link_pm(left_node, right_invnode))
                    graph_file.write(format_gfa_link_mp(left_invnode, right_node))

                    # Associate alternative links to sv_id (INV)
                    link_id_1 = (left_node, "+", right_invnode, "-")
                    link_id_2 = (left_invnode, "-", right_node, "+")

                    for link_id in [link_id_1, link_id_2]:

                        link_key = get_link_key(link_id)

                        if link_key not in d_link_sv.keys():
                            d_link_sv[link_key] = []
                        
                        d_link_sv[link_key].append((":".join([chrom, sv_id]), 1))

            elif sv_type == "BND":
                
                left_coords, right_coords = parse_BND_id(chrom, sv_id)

                ## DEBUG
                # if 0 in [int(left_coords[1]), int(right_coords[1])]:
                #     print(sv_id, left_coords, right_coords)
                ##

                if left_coords[2] == "-":
                    left_node = find_node_by_start(left_coords[0], int(left_coords[1]), all_nodes[left_coords[0]])
                else:
                    left_node = find_node_by_end(left_coords[0], int(left_coords[1]), all_nodes[left_coords[0]])

                if right_coords[2] == "+":
                    right_node = find_node_by_start(right_coords[0], int(right_coords[1]), all_nodes[right_coords[0]])
                else:
                    right_node = find_node_by_end(right_coords[0], int(right_coords[1]), all_nodes[right_coords[0]])
                
                if not None in [left_node, right_node]:

                    # Create alternative link and associate with sv_id
                    if left_coords[2] == "-":
                        graph_file.write(format_gfa_link_mp(left_node, right_node))
                        link_id = (left_node, "-", right_node, "+")

                    elif right_coords[2] == "-":
                        graph_file.write(format_gfa_link_pm(left_node, right_node))
                        link_id = (left_node, "+", right_node, "-")

                    else:
                        graph_file.write(format_gfa_link_pp(left_node, right_node, 1))
                        link_id = (left_node, "+", right_node, "+")
                    
                    #Associate with sv_id
                    link_key = get_link_key(link_id)

                    if link_key not in d_link_sv.keys():
                        d_link_sv[link_key] = []

                    d_link_sv[link_key].append((":".join([chrom, sv_id]), 1))
                
                else:
                    print(f"Warning: no alternative link defined for {sv_id} (ID: {d_sv_ID[sv_id]})")

    graph_file.close()


    # Output d_link_sv as json file
    with open("svs_edges.json", 'w') as file:
        file.write(json.dumps(d_link_sv, sort_keys=True, indent=4))

def get_link_key(link_id):
    return "@".join(link_id)

def format_gfa_link_pp(node1, node2, rank):
    # print("\t".join(["L", node1, "+", node2, "+", "0M", "SR:i:"+str(rank)]) + "\n")
    return "\t".join(["L", node1, "+", node2, "+", "0M"]) + "\n"
def format_gfa_link_pm(node1, node2):
    rank = 1
    return "\t".join(["L", node1, "+", node2, "-", "0M"]) + "\n"
def format_gfa_link_mp(node1, node2):
    rank = 1
    return "\t".join(["L", node1, "-", node2, "+", "0M"]) + "\n"
def format_gfa_link_mm(node1, node2):
    rank = 1
    return "\t".join(["L", node1, "-", node2, "-", "0M"]) + "\n"

def format_gfa_node(node_id, seq):
    '''Convert node info in GFA format. For node_id format pattern, see function format_node_id().'''
    return "\t".join(["S", node_id, seq]) + "\n"
def format_gfa_path(name, nodes, lengths):
    return "\t".join(["P", name, nodes, lengths]) + "\n"

def format_node_id(chrom, coord1, coord2):
    # Add +1 to start and end positions to use 1-indexed positions (compatible with 1-indexed positions of VCF file format)
    return ":".join([ str(chrom), "-".join([ str(coord1+1), str(coord2+1) ]) ])
def format_altnode_id(chrom, coord1, count):
    return ":".join([ str(chrom), str(coord1)+f".{str(count)}" ])

def find_node_by_start(Tchr, pos, node_list):
    for node in node_list:
        chrom, coords = node.split(":")
        start_pos = int(coords.split("-")[0])
        if chrom == Tchr and start_pos == pos:
            return node
    #Test
    # sys.exit("Node not found, incorrect start position: " + str(pos))
    print(f"Warning: looked for nonexistant node starting at {str(pos)} on {Tchr}")
    return None

def find_node_by_end(Tchr, pos, node_list):
    for node in node_list:
        chrom, coords = node.split(":")
        end_pos = int(coords.split("-")[1])
        if chrom == Tchr and end_pos == pos:
            return node
    #Test
    # sys.exit("Node not found, incorrect end position: " + str(pos))
    print(f"Warning: looked for nonexistant node ending at {str(pos)} on {Tchr}")
    return None

def format_DEL_id(pos, end):
    return '-'.join(["DEL", str(pos), str(end)])

def format_INS_id(pos, ins_count):
    return '-'.join(["INS", str(pos), str(ins_count)])

def format_INV_id(pos, end):
    return '-'.join(["INV", str(pos), str(end)])

def format_BND_id(pos, alt):
    '''
    There are 4 types of BND possible (https://samtools.github.io/hts-specs/VCFv4.2.pdf):
    
    - t[p[: [+++t] -> [p+++]
    - t]p]: [+++t] -> [---p] (reverse comp piece extending left of p is joined after t)
    - ]p]t: [+++p] -> [t+++]
    - [p[t: [p---] -> [t+++] (reverse comp piece extending right of p is joined before t)

    where t is the nucleotide indicated in the REF field (at POS in CHROM),
    and p is the ending position of bkpt (chrom2:pos2)
    '''

    # BND type: t[p[ or [p[t
    if "[" in alt:
        parsed_alt = list(filter(bool, alt.split("[")))

        # t[p[: piece extending to the right of p is joined after t
        if ":" in parsed_alt[1]:
            t = parsed_alt[0]

        # [p[t: reverse comp piece extending right of p is joined before t
        else:
            t = parsed_alt[1]

    # BND type: t]p] or ]p]t
    elif "]" in alt:
        parsed_alt = list(filter(bool, alt.split("]")))

        # t]p]: reverse comp piece extending left of p is joined after t
        if ":" in parsed_alt[1]:
            t = parsed_alt[0]
        
        # ]p]t: piece extending to the left of p is joined before t
        else:
            t = parsed_alt[1]

    else:
        return "BND-format"

    # if len(t) > 1:
    #     return "BND-format"

    alt = alt.replace(t, pos)

    return '-'.join(["BND", alt])

def parse_BND_id(chrom, bnd_id):

    right_chrom = None
    alt = bnd_id.split("BND-")[1]

    if "[" in alt:
        
        alt = list(filter(bool, alt.split("[")))

        # t[p[: piece extending to the right of p is joined after t
        if ":" in alt[1]:

            t = alt[0]
            p = alt[1]

            left_chrom = chrom
            left_pos = t
            left_strand = "+"

            right_chrom = p.split(":")[0]
            right_pos = p.split(":")[1]
            right_strand = "+"

        # [p[t: reverse comp piece extending right of p is joined before t
        elif ":" in alt[0]:
            
            t = alt[1]
            p = alt[0]

            left_chrom = p.split(":")[0]
            left_pos = p.split(":")[1]
            left_strand = "-"

            right_chrom = chrom
            right_pos = t
            right_strand = "+"

    elif "]" in alt:

        alt = list(filter(bool, alt.split("]")))

        # t]p]: reverse comp piece extending left of p is joined after t
        if ":" in alt[1]:

            t = alt[0]
            p = alt[1]
            
            left_chrom = chrom
            left_pos = t
            left_strand = "+"

            right_chrom = p.split(":")[0]
            right_pos = p.split(":")[1]
            right_strand = "-"

        # ]p]t: piece extending to the left of p is joined before t
        elif ":" in alt[0]:

            t = alt[1]
            p = alt[0]
            
            left_chrom = p.split(":")[0]
            left_pos = p.split(":")[1]
            left_strand = "+"

            right_chrom = chrom
            right_pos = t
            right_strand = "+"

    if right_chrom is None:
        return None, None
    else:
        return [left_chrom, int(left_pos), left_strand], [right_chrom, int(right_pos), right_strand]


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
        inVCF, inFA, outGFA = parse_arguments(sys.argv[1:])
        construct_gfa(inVCF, inFA, outGFA)