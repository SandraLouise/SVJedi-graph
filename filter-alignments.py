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
        "-O", 
        "--dover", metavar="<min_breakpoint_overlap>", nargs=1, 
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
    d_over = args.dover

    #=========================================================
    # 1. Load graph info
    #=========================================================
    # 1.1. Edges - SVs associations
    #---------------------------------------------------------
    with open("svs_edges.json", "r") as link_sv_json:
        d_link_sv = json.load(link_sv_json)
        # key = link_id = "left_node_id@strand@right_node_id@strand"
        # value = list((sv_id, allele))

    #---------------------------------------------------------
    # 1.2. Length of alt nodes
    #---------------------------------------------------------
    alt_node_len = {}
    
    with open(gfa_file, "r") as graph_file:
        for line in graph_file:

            #Process alternative nodes (INS nodes)
            if line.startswith("S") and "." in line.split("\t")[1].split(":")[-1]:

                node = line.split("\t")[1]
                # Save node length
                alt_node_len[node] = len(line.rstrip().split("\t")[2])

    #=========================================================
    # 2. Read alignments and fill dictionary of informative alignments
    #=========================================================

    dict_of_informative_aln = dict()
    list_of_identities = list()


    with open(gaf_file) as aln_file:
        for line in aln_file:

            aln = read_gaf_line(line.rstrip())
            identity_saved = False
    
            # Process path 
            target_nodes = extract_nodes(aln["Tid"])

            # Filter < 2 nodes alns (no breakpoint overlap)
            if len(target_nodes) < 2:
                continue
            
            target_links = get_aln_links(aln["Tid"], target_nodes)

            # Get overlapped SVs for each link in the alignment
            for link in target_links:

                link_key = get_link_key(link)
                rev_link_key = get_link_key(reverse_link(link))

                # Check for presence of both link and rev_link in d_link_sv
                key_to_check = []
                for key in [link_key, rev_link_key]:
                    if key in d_link_sv.keys():
                        key_to_check.append(key)

                for key in key_to_check:

                    # Check each overlapped SV
                    for (sv_id, allele) in d_link_sv[key]:
                        
                        # BREAKPOINT FILTER
                        if check_bkpt_overlap(link, target_nodes, d_over, aln, alt_node_len):
                            
                            # Add aln to dict
                            # here sv_id = {chrom}:{sv_type}-{pos description for this sv type}
                            sv_type = sv_id.split(":")[1].split("-")[0]
                            dict_sv_id = format_dict_sv_id(sv_id, sv_type)

                            if dict_sv_id not in dict_of_informative_aln.keys():
                                dict_of_informative_aln[dict_sv_id] = [[], []]

                            dict_of_informative_aln[dict_sv_id][allele].append(line.split("cg:Z:")[0])


    #Stats
    aln_nb = 0
    for sv in dict_of_informative_aln:
        aln_nb += len(dict_of_informative_aln[sv][0]) + len(dict_of_informative_aln[sv][1])

    with open(output_aln_dict, 'w') as file:
        file.write(json.dumps(dict_of_informative_aln, sort_keys=True, indent=4))
    # print(str(aln_nb), "informative sv-based aln saved")
    # print(len(removed_alns), " alns removed")

    # print(dict_of_informative_aln)

def get_link_key(link_id):
    return "@".join(link_id)

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

def get_aln_links(path, nodes):

    # Get nodes strands
    strands = [None]*len(nodes)
    for i in range(len(nodes)):

        if path.split(nodes[i])[0][-1] == ">":
            strands[i] = "+"
        else:
            strands[i] = "-"

    # Get aln links
    links = [None]*(len(nodes)-1)

    for i in range(1, len(nodes)):
        link_id = (nodes[i-1], strands[i-1], nodes[i], strands[i])

        links[i-1] = link_id

    return links

def reverse_link(link):

    rev_strand = { "+" : "-", "-" : "+" }

    return (link[2], rev_strand[link[3]], link[0], rev_strand[link[1]])

def invert_orient(orient):
    d = {">" : "<", "<" : ">"}
    return d[orient]

def format_DEL_INS_id(sv_type, pos, end):
    return '-'.join([sv_type, str(pos), str(end)])
def format_INV_id(pos, end):
    return '-'.join(["INV", str(pos), str(end)])
def format_BND_id(pos, end_chr, end_pos):
    return '-'.join(["BND", str(pos), ":".join([end_chr, end_pos])])

def format_dict_sv_id(aln_dict_sv_id, sv_type):
    '''
    Became useless after modification of predict-genotype.py to 
    match sv_id in dict_of_informative_aln.
    Kept to keep track of sv_id format.

    Format in dict_of_informative_aln: 
                {chrom} : {sv_id as formatted by functions of contruct-graph.py}
    - DEL:      {chrom} : DEL - {pos} - {end}
    - INS:      {chrom} : INS - {pos} - {ins_count} 
    - INV:      {chrom} : INV - {pos} - {end}
    - BND:      {chrom} : BND - {alt}
                with alt = {pos}[{chrom2}:{pos2[
                     or    {pos}]{chrom2}:{pos2]
                     or    [{chrom2}:{pos2[{pos}
                           ]{chrom2}:{pos2]{pos}
    '''

    return aln_dict_sv_id

def check_bkpt_overlap(link, aln_nodes, d_over, aln, alt_node_len):

    unaligned_start = aln["Ts"]
    unaligned_end = aln["Tlen"] - aln["Te"] - 1

    bkpt_leftNode = link[0]
    bkpt_rightNode = link[2]

    left_overlap = False
    right_overlap = False

    left_overlap = sum(get_node_len(node, alt_node_len) for node in aln_nodes[:aln_nodes.index(bkpt_leftNode)+1]) - unaligned_start >= d_over

    right_overlap = sum(get_node_len(node, alt_node_len) for node in aln_nodes[aln_nodes.index(bkpt_rightNode):]) - unaligned_end >= d_over

    return left_overlap and right_overlap

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
    #if alternative node
    if "." in node_id.split(":")[-1]:
        return node_id.split(":")[-1].split(".")[0]
    #if reference node
    else:
        return node_id.split(":")[-1].split("-")[0]
def get_node_end(node_id):
    #if alternative node
    if "." in node_id.split(":")[-1]:
        return node_id.split(":")[-1].split(".")[0]
    #if reference node
    else:
        return node_id.split(":")[-1].split("-")[1]
def get_node_len(node_id, alt_node_len):
    #if alternative node
    if "." in node_id.split(":")[-1]:
        return alt_node_len[node_id]
    #if reference node
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
