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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""

import sys
import argparse
import re
import json

def main(args):

    parser = argparse.ArgumentParser(
        description="---"
    )    
    parser.add_argument(
        "-a", 
        "--gaf", metavar="<align_file>", nargs=1, 
        help="align file in gaf format", 
        required=True
    )
    parser.add_argument(
        "-g", 
        "--gfa", metavar="<graph_file>", nargs=1, 
        help="variant graph in gfa format", 
        required=True
    )

    parser.add_argument(
        "-i", 
        "--gfainfo", metavar="<gfa_info>", nargs=1, 
        help="gfa info", 
        required=False
    )

    parser.add_argument(
        "-L", 
        "--ladj", metavar="<len_adjacent_seq>", nargs=1, 
        required=False, 
        default=5000
    )
    parser.add_argument(
        "-O", 
        "--dover", metavar="<min_breakpoint_overlap>", nargs=1, 
        required=False, 
        default=100
    )
    parser.add_argument(
        "-E", 
        "--dend", metavar="<max_dist_semiglobal>", nargs=1, 
        required=False, 
        default=100
    )

    parser.add_argument(
        "-o", 
        "--outputDir", 
        metavar="<outputDirectory>", 
        type=str,
        required=False)

    parser.add_argument(
        "-p", 
        "--prefix", 
        metavar="<prefix", 
        type=str,
        required=False)

    args = parser.parse_args()
        
    if args.prefix:
        output_aln_dict = args.prefix + '_informative_aln' + '.json'
        gfa_info_dict = args.prefix + "_gfa_info.json"
    else:
        output_aln_dict = 'informative_aln.json'
        gfa_info_dict = "gfa_info.json"
    
    if args.outputDir:
        output_aln_dict = "/".join([args.outputDir, output_aln_dict])

    gfa_file = args.gfa[0]
    gaf_file = args.gaf[0]
    l_adj = args.ladj
    d_over = args.dover
    d_end = args.dend


    if args.gfainfo is None:
        #1. Listing SVs to genotype
        # dict_SVs = {} # d = {chrom : [sv_id, ...]}
        region_svs = {}
        which_chrom = {}
        start_end_chrom = {}
        rename_node = {}

        with open(gfa_file, "r") as graph_file:
            for line in graph_file:
                if line.startswith("#"):
                    region_svs[line[1:].split("\t")[0]] = list(line.rstrip().split("\t")[1].split(";"))

                elif line.startswith("P"):
                    region, nodes, __ = line.split("\t")[1:]
                    # chrom = region

                    # Add nodes and region to which_chrom
                    nodes = extract_nodes(nodes)
                    for node in nodes:
                        which_chrom[node] = region

                    # Add first and last node to start_end_chrom
                    start_end_chrom[region] = [nodes[0], nodes[-1]]

                elif line.startswith("S"):
                    node = line.split("\t")[1]
                    if node[-1] == "a":
                        for n_id in list(which_chrom.keys()):
                            if n_id.startswith(node[:-1]):
                                which_chrom[node] = which_chrom[n_id]
                                break
                            elif n_id.startswith("{}:{}".format(node.split(":")[0], str(int(node.split(":")[1][:-1]) - 1))):
                                new_id = "{}:{}a".format(node.split(":")[0], str(int(node.split(":")[1][:-1]) - 1))
                                which_chrom[new_id] = which_chrom[n_id]
                                rename_node[node] = new_id
                                # print(node, n_id, new_id)
                                break
                ###########################
        
        with open(gfa_info_dict, 'w') as file:
            gfa_info = {"which_chrom" : which_chrom, "rename_node" : rename_node, "start_end_chrom" : start_end_chrom}
            file.write(json.dumps(gfa_info, sort_keys=True, indent=4))
    
    else:
        with open(args.gfainfo[0], "r") as file:
            gfa_info = json.load(file)
            which_chrom = gfa_info["which_chrom"]
            start_end_chrom = gfa_info["start_end_chrom"]
            rename_node = gfa_info["rename_node"]


    #2. Reading alignments and filling dictionary of informative alignments

    dict_of_informative_aln = dict()
    list_of_identities = list()
    # dict = {sv_id : [[aln on ref allele], [aln on alt allele]]}


    with open(gaf_file) as aln_file:
        for line in aln_file:

            aln = read_gaf_line(line.rstrip())
            identity_saved = False
    
            # Process path 
            target_nodes = extract_nodes(aln["Tid"])

            #filter < 2 nodes alns (no breakpoint overlap)
            if len(target_nodes) < 2:
                continue
            
            # SEMI-GLOBALITY FILTER
            multi_chrom = False

            target_chrom = target_nodes[0].split(":")[0]
            targets = [target_chrom]
                
            for node in target_nodes[1:]:

                if node.split(":")[0] != target_chrom:
                    multi_chrom = True
                    targets.append(which_chrom[node])
                    
            # Inter-chromosomal BNDs
            if multi_chrom:
                target_chrom = targets
                #TODO: handle multi target-regions (needed for TRANSLOCATIONS)

            # Convert path_based aln to sv-based aln
            target_svs = get_target_svs(target_nodes, target_chrom, multi_chrom, region_svs)

            for sv_id in target_svs: #sv_id = sv_type + "-" + pos + "-" + end

                sv_type = sv_id.split(":")[1].split("-")[0]

                breakpoints = get_breakpoints(sv_id, sv_type)

                #identify supported allele
                allele = get_allele(breakpoints, aln["Tid"], target_nodes)

                if allele is None:
                    continue

                # BREAKPOINT FILTER
                overlapped_bkpts = overlap_breakpoints(breakpoints, sv_id, allele, target_nodes, d_over, aln)

                # Add aln to dict and save identity
                dict_sv_id = format_dict_sv_id(sv_id, sv_type)

                if dict_sv_id not in dict_of_informative_aln.keys():
                    dict_of_informative_aln[dict_sv_id] = [[], []]

                for bkpt in overlapped_bkpts:
                    dict_of_informative_aln[dict_sv_id][allele].append(line.split("cg:Z:")[0])

                if not identity_saved:
                    list_of_identities.append(aln["Aid"])
                    identity_saved = True


    with open(output_aln_dict, 'w') as file:
        file.write(json.dumps(dict_of_informative_aln, sort_keys=True, indent=4))



def read_gaf_line(line):
    Qid, Qlen, Qs, Qe = line.split("\t")[:4]
    Tid, Tlen, Ts, Te = line.split("\t")[5:9]
    Am, Alen, Aq = line.split("\t")[9:12]
    
    d = {"Qid" : Qid, "Qlen" : int(Qlen), "Qs" : int(Qs), "Qe" : int(Qe),
         "Tid" : Tid, "Tlen" : int(Tlen), "Ts" : int(Ts), "Te" : int(Te), 
         "Am" : int(Am), "Alen" : int(Alen), "Aq" : int(Aq)}
    
    if "id:f:" in line:
        d["Aid"] = float(line.split("id:f:")[-1].split("\t")[0])
    else:
        d["Aid"] = d["Am"] / d["Alen"]
        
    return d

def invert_orient(orient):
    d = {">" : "<", "<" : ">"}
    return d[orient]

def get_allele(bkpts, path, nodes):

    for node in nodes:
        if node.endswith("a"):
            return 1

    for allele in [0,1]:
        for bkpt in bkpts[allele]:
            # print(bkpt)
        
            for i in range(len(nodes)-1):
           
                # bkpt = [chrom:pos_left, chrom:pos_right]
                chr1, pos1, orient1 = bkpt[0].split(":")
                chr2, pos2, orient2 = bkpt[1].split(":")

                # Check if bkpt link in aln path
                ## Left node
                if orient1 == ">":
                    left = all([nodes[i].startswith(chr1), get_node_end(nodes[i]) == pos1, path.split(nodes[i])[0][-1] == orient1])
                else:
                    left = all([nodes[i].startswith(chr1), get_node_start(nodes[i]) == pos1, path.split(nodes[i])[0][-1] == orient1])

                ## Right node
                if orient2 == ">":
                    right = all([nodes[i+1].startswith(chr2), get_node_start(nodes[i+1]) == pos2, path.split(nodes[i+1])[0][-1] == orient2])
                else:
                    right = all([nodes[i+1].startswith(chr2), get_node_end(nodes[i+1]) == pos2, path.split(nodes[i+1])[0][-1] == orient2])

                # Check in bkpt link in aln path (reverse)
                ## Left node
                if orient1 == ">":
                    left_rev = all([nodes[i+1].startswith(chr1), get_node_end(nodes[i+1]) == pos1, path.split(nodes[i+1])[0][-1] == invert_orient(orient1)])
                else:
                    left_rev = all([nodes[i+1].startswith(chr1), get_node_start(nodes[i+1]) == pos1, path.split(nodes[i+1])[0][-1] == invert_orient(orient1)])
                
                ## Right node
                if orient2 == ">":
                    right_rev = all([nodes[i].startswith(chr2), get_node_start(nodes[i]) == pos2, path.split(nodes[i])[0][-1] == invert_orient(orient2)])
                else:
                    right_rev = all([nodes[i].startswith(chr2), get_node_end(nodes[i]) == pos2, path.split(nodes[i])[0][-1] == invert_orient(orient2)])

                if (left and right) or (left_rev and right_rev):
                    return allele

    
    return None


def format_DEL_INS_id(sv_type, pos, end):
    return '-'.join([sv_type, str(pos), str(end)])
def format_INV_id(pos, end):
    return '-'.join(["INV", str(pos), str(end)])
def format_BND_id(pos, end_chr, end_pos):
    return '-'.join(["BND", str(pos), ":".join([end_chr, end_pos])])

def format_dict_sv_id(graph_sv_id, sv_type):
    '''
    Graph sv_id:
    - DEL: chrom + "_" + sv_type + "-" + pos + "-" + end
    - INS: chrom + "_" + sv_type + "-" + pos + "-" + end
    - INV: chrom + "_" + sv_type + "-" + pos + "-" + end
    - BND: chrom + "_" + sv_type + "-" + pos + "-" + end_chr + ":" + end_pos

    Dict sv_id:
    - DEL/INS/INV: chrom + "_" + pos + "-" + end
    - BND: ?
    '''
    # print(graph_sv_id)

    return graph_sv_id.replace(sv_type + "-", "")

def get_breakpoints(sv_id, svtype):
    '''Returns all breakpoints coordinates (as a list of one duple per bkpt) for a SV (sorted by allele)'''
    #Considering 2 positions possible for left coord (SV start) because of particularities of real SV evaluation set

    bkpts = [[], []] # [[allele 0], [allele 1]]
    forward = ">"
    reverse = "<"
    
    if svtype != "BND":

        chrom, info = sv_id.split(":")
        __, pos, end = info.split("-")

        if svtype == "DEL":
            bkpts[1] = [(":".join([chrom, pos, forward]), ":".join([chrom, str(int(end)+1), forward]))]
            bkpts[0] = [(":".join([chrom, pos, forward]), ":".join([chrom, str(int(pos)+1), forward])), 
                        (":".join([chrom, end, forward]), ":".join([chrom, str(int(end)+1), forward]))]
        
        elif svtype == "INS":
            end = str(int(pos) + 1)
            bkpts[0] = [(":".join([chrom, pos, forward]), ":".join([chrom, end, forward]))]
            bkpts[1] = [(":".join([chrom, pos, forward]), ":".join([chrom, end+"a", forward])), 
                        (":".join([chrom, end+"a", forward]), ":".join([chrom, end, forward]))]
        
        elif svtype == "INV":
            bkpts[0] = [(":".join([chrom, pos, forward]), ":".join([chrom, str(int(pos)+1), forward])), 
                        (":".join([chrom, end, forward]), ":".join([chrom, str(int(end)+1), forward]))]
            bkpts[1] = [(":".join([chrom, pos, forward]), ":".join([chrom, end, reverse])), 
                        (":".join([chrom, str(int(pos)+1), reverse]), ":".join([chrom, str(int(end)+1), forward]))]

    else:
        chrom, bnd_id_p1, bnd_id_p2 = sv_id.split(":")
        left_coords, right_coords = parse_BND_id(chrom, ":".join([bnd_id_p1, bnd_id_p2]))

        # Reference breakpoints
        if left_coords[2] == forward:
            bkpts[0].append((":".join(left_coords), ":".join([left_coords[0], str(int(left_coords[1])+1), forward])))
        else:
            bkpts[0].append((":".join([left_coords[0], str(int(left_coords[1])-1), forward]), ":".join(left_coords)))
        
        if right_coords[2] == forward:
            bkpts[0].append((":".join([right_coords[0], str(int(right_coords[1])-1), forward]), ":".join(right_coords)))
        else:
            bkpts[0].append((":".join(right_coords), ":".join([right_coords[0], str(int(right_coords[1])+1), forward])))
        
        # Alternative breakpoint
        bkpts[1].append((":".join(left_coords), ":".join(right_coords)))

    return bkpts

def overlap_breakpoints(bkpts, sv_id, allele, aln_nodes, d_over, aln):

    over_bkpt = []

    for bkpt in bkpts[allele]:
        if check_single_breakpoint(sv_id, bkpt, aln_nodes, d_over, aln):
            over_bkpt.append(bkpt)

    return over_bkpt

def check_single_breakpoint(sv_id, bkpt, aln_nodes, d_over, aln):
    bkpt_leftChrom, bkpt_leftCoord, bkpt_leftStrand = bkpt[0].split(":")
    bkpt_rightChrom, bkpt_rightCoord, bkpt_rightStrand = bkpt[1].split(":")
    unaligned_start = aln["Ts"]
    unaligned_end = aln["Tlen"] - aln["Te"] - 1

    #Check presence of breakpoint nodes
    found_leftNode = False
    found_rightNode = False
    for node in aln_nodes:

        # Look for left node of bkpt
        if bkpt_leftStrand == ">":
            if get_node_end(node) == bkpt_leftCoord and get_node_chrom(node) == bkpt_leftChrom:
                found_leftNode = True
                bkpt_leftNode = node
        else:
            if get_node_start(node) == bkpt_leftCoord and get_node_chrom(node) == bkpt_leftChrom:
                found_leftNode = True
                bkpt_leftNode = node
        
        # Look for right node of bkpt
        if bkpt_rightStrand == ">":
            if get_node_start(node) == bkpt_rightCoord and get_node_chrom(node) == bkpt_rightChrom:
                found_rightNode = True
                bkpt_rightNode = node
        else:
            if get_node_end(node) == bkpt_rightCoord and get_node_chrom(node) == bkpt_rightChrom:
                found_rightNode = True
                bkpt_rightNode = node

    if not (found_leftNode and found_rightNode):
        return False

    #Check breakpoint overlap
    left_overlap = False
    right_overlap = False

    #Forward alignment
    if aln_nodes.index(bkpt_leftNode) < aln_nodes.index(bkpt_rightNode):
        left_overlap = sum(get_node_len(node, sv_id) for node in aln_nodes[:aln_nodes.index(bkpt_leftNode)+1]) - unaligned_start >= d_over
        right_overlap = sum(get_node_len(node, sv_id) for node in aln_nodes[aln_nodes.index(bkpt_rightNode):]) - unaligned_end >= d_over

    
    #Reverse alignment
    else:
        left_overlap = sum(get_node_len(node, sv_id) for node in aln_nodes[aln_nodes.index(bkpt_leftNode):]) - unaligned_end >= d_over
        right_overlap = sum(get_node_len(node, sv_id) for node in aln_nodes[:aln_nodes.index(bkpt_rightNode)+1]) - unaligned_start >= d_over

    return left_overlap and right_overlap   

def is_semiglobal_aln(aln, target_nodes, start_end_region, d_end):
    # a_coords = [[target_nodes[0], int(p_start)], [target_nodes[-1], int(p_len)-int(p_end)-1]]
    # a_coords_for_semiglob = [[a_start_node, unaligned_left], [a_end_node, unaligned_right]]
    # list_component_first_last_node = [component first node, component last node]

    glob_region_left_f = all([start_end_region[0] == target_nodes[0], aln["Ts"] <= d_end])
    glob_region_left_r = all([start_end_region[0] == target_nodes[-1], aln["Tlen"] - aln["Te"] - 1 <= d_end])
    glob_region_left = any([glob_region_left_f, glob_region_left_r])

    glob_region_right_f = all([start_end_region[1] == target_nodes[-1], aln["Tlen"] - aln["Te"] - 1 <= d_end])
    glob_region_right_r = all([start_end_region[1] == target_nodes[0], aln["Ts"] <= d_end <= d_end])
    glob_region_right = any([glob_region_right_f, glob_region_right_r])

    return any([(glob_region_left) and (glob_region_right), #allele inside read
                (aln["Qs"] <= d_end) and (aln["Qlen"] - d_end <= aln["Qe"]), #read inside allele
                (aln["Qs"] <= d_end) and (glob_region_right), #read left aligned
                (glob_region_left) and (aln["Qlen"] - d_end <= aln["Qe"])]) #read right aligned

def get_target_svs(target_nodes, target_chrom, multi_chrom, region_svs):

    target_svs = []

    def get_left_right_ends(node_list):
        left_node_start = get_node_start(node_list[0])
        if left_node_start.endswith("a"):
            left_node_start = int(left_node_start[:-1])
        else:
            left_node_start = int(left_node_start)

        right_node_end = get_node_end(node_list[-1])
        if right_node_end.endswith("a"):
            right_node_end = int(right_node_end[:-1])
        else:
            right_node_end = int(right_node_end)
        
        return left_node_start, right_node_end


    if not multi_chrom:

        if target_chrom not in region_svs.keys():
            return target_svs

        left_node_start, right_node_end = get_left_right_ends(target_nodes)

        for sv_id in region_svs[target_chrom]:

            if sv_id.startswith("BND-"):
                left_coords, right_coords = parse_BND_id(target_chrom, sv_id)
                sv_start = left_coords[1]
                sv_end = right_coords[1]

            else:
                __, sv_start, sv_end = sv_id.split("-")

            if any([min(left_node_start, right_node_end) <= int(sv_start) <= max(left_node_start, right_node_end), min(left_node_start, right_node_end) <= int(sv_end) <= max(left_node_start, right_node_end)]):
                target_svs.append(":".join([target_chrom, sv_id])) #add chrom to id for later save in aln dict

    else:

        for chrom in target_chrom:
            if chrom not in region_svs.keys():
                return target_svs

        for chrom in target_chrom:

            chrom_nodes = [node for node in target_nodes if node.startswith(chrom)]
            left_node_start, right_node_end = get_left_right_ends(chrom_nodes)

            for sv_id in region_svs[chrom]:
                if sv_id.startswith("BND-"):
                    left_coords, right_coords = parse_BND_id(target_chrom, sv_id)

                    if left_coords[0] == chrom:
                        sv_pos = left_coords[1]
                        chrom_bis, sv_pos_bis = right_coords[:2]
                        chrom_bis_nodes = [node for node in target_nodes if node.startswith(chrom_bis)]
                        left_node_start_bis, right_node_end_bis = get_left_right_ends(chrom_bis_nodes)
                    
                    else:
                        sv_pos = right_coords[1]
                        chrom_bis, sv_pos_bis = left_coords[:2]
                        chrom_bis_nodes = [node for node in target_nodes if node.startswith(chrom_bis)]
                        left_node_start_bis, right_node_end_bis = get_left_right_ends(chrom_bis_nodes)
                    
                    if any([min(left_node_start, right_node_end) <= int(sv_pos) <= max(left_node_start, right_node_end), min(left_node_start_bis, right_node_end_bis) <= int(sv_pos_bis) <= max(left_node_start_bis, right_node_end_bis)]):
                        target_svs.append(":".join([chrom, sv_id]))
                
                else:
                    __, sv_start, sv_end = sv_id.split("-")

                    if any([min(left_node_start, right_node_end) <= int(sv_start) <= max(left_node_start, right_node_end), min(left_node_start, right_node_end) <= int(sv_end) <= max(left_node_start, right_node_end)]):
                        target_svs.append(":".join([chrom, sv_id]))

        #TODO: handle multi-chrom BNDs

    return target_svs

def parse_BND_id(chrom, bnd_id):

    right_chrom = None
    alt = bnd_id.split("BND-")[1]

    if "[" in alt:
        alt = list(filter(bool, alt.split("[")))
        if ":" in alt[1]:
            # t[p[ : piece extending to the right of p is joined after t
            right_chrom = alt[1].split(":")[0]
            right_pos = alt[1].split(":")[1]
            right_strand = ">"
            left_chrom = chrom
            left_pos = alt[0]
            left_strand = ">"

        elif ":" in alt[0]:
            # [p[t : reverse comp piece extending right of p is joined before t
            left_chrom = alt[0].split(":")[0]
            left_pos = alt[0].split(":")[1]
            left_strand = "<"
            right_chrom = chrom
            right_pos = alt[1]
            right_strand = ">"

    elif "]" in alt:
        alt = list(filter(bool, alt.split("]")))
        if ":" in alt[1]:
            # t]p] : reverse comp piece extending left of p is joined after t
            right_chrom = alt[1].split(":")[0]
            right_pos = alt[1].split(":")[1]
            right_strand = "<"
            left_chrom = chrom
            left_pos = alt[0]
            left_strand = ">"

        elif ":" in alt[0]:
            # ]p]t piece extending to the left of p is joined before t
            left_chrom = alt[0].split(":")[0]
            left_pos = alt[0].split(":")[1]
            left_strand = ">"
            right_chrom = chrom
            right_pos = alt[1]
            right_strand = ">"

    if right_chrom is None:
        return None, None
    else:
        return (left_chrom, left_pos, left_strand), (right_chrom, right_pos, right_strand)

def get_node_chrom(node_id):
    return node_id.split(":")[0]

def get_node_start(node_id):
    #node_id = chrom : start - end
    if node_id.endswith("a"):
        return node_id.split(":")[1]
    else:
        return node_id.split(":")[1].split("-")[0]
def get_node_end(node_id):
    if node_id.endswith("a"):
        end = str(int(get_node_start(node_id)[:-1]) + 1) + "a"
    else:
        end = node_id.split(":")[1].split("-")[1]
    return end
def get_node_len(node_id, sv_id):
    if node_id.endswith("a"):
        return int(sv_id.split("-")[-1])
    else:
        return int(get_node_end(node_id)) - int(get_node_start(node_id)) + 1

def extract_nodes(p, *args, **kwargs):
    """
    Return nodes in list format from nodes in str format. 
    Can take path orientation into account if input is gaf path.

    Args:
        str_node_list (str) : node list in str, with nodes separated by special characters
        orient (bool) : checks orientation of nodes in aln path (ONLY WITH GAF PATH)

    Ex: '1+,2+,3+,4+' => [1, 2, 3, 4]   (path in gfa)
           '>1>2>3>4' => [1, 2, 3, 4]   (path in gaf)
    """
    # check_orient = kwargs.get('orient', False)

    #Path in GAF
    if p[0] in [">", "<"]:
        nodes = [s for s in re.split(r'[<>]', p) if s]

    #Path in GFA
    else:
        nodes = [s[:-1] for s in re.split(r'[,]', p) if s]

    return nodes


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])
