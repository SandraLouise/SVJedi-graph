#!/usr/bin/env python3

import sys
import argparse
import re
import json
import statistics as stats

def main(args):

    parser = argparse.ArgumentParser(
        description="---"
    )    
    parser.add_argument(
        "-a", "--gaf", metavar="<align_file>", nargs=1, help="align file in gaf format", required=True
    )
    parser.add_argument(
        "-g", "--gfa", metavar="<graph_file>", nargs=1, help="variant graph in gfa format", required=True
    )

    parser.add_argument(
        "-n", "--sv_nodes", metavar="<sv_node_dict>", nargs=1, help="sv nodes info", required=False
    )

    parser.add_argument(
        "-L", "--ladj", metavar="<len_adjacent_seq>", nargs=1, required=False, default=5000
    )
    parser.add_argument(
        "-O", "--dover", metavar="<min_breakpoint_overlap>", nargs=1, required=False, default=100
    )
    parser.add_argument(
        "-E", "--dend", metavar="<max_dist_semiglobal>", nargs=1, required=False, default=100
    )

    parser.add_argument(
        "-v", 
        "--version", 
        metavar="<versionID", 
        type=str,
        required=False)

    args = parser.parse_args()
    
    if args.sv_nodes is not None:
        pass #load dict_of_DEL_sv_nodes from json

    if args.version:
        output_aln_dict = 'informative_aln_' + args.version + '.json'
    else:
        output_aln_dict = 'informative_aln.json'

    gfa_file = args.gfa[0]
    gaf_file = args.gaf[0]
    l_adj = args.ladj
    d_over = args.dover
    d_end = args.dend


    #1. Preparing dict of SV nodes
    print("Reading gfa file...")

    if args.sv_nodes is not None:
        pass #don't recreate dict (real time save or not?)

    dict_of_sv_nodes = dict()
    #dict = {sv_id : [[sv region nodes], [sv nodes]]}
    dict_of_nodes_len = dict()
    #dict = {node_id : node_len}
    
    #1.1. Writing node list for each SV 

    with open(gfa_file) as graph_file:
        for line in graph_file:  

            if line[0] == "S":
                #Node line
                __, node_id, node_seq = line.rstrip().split("\t")
                dict_of_nodes_len[int(node_id)] = len(node_seq)

            if line[0] == "P":
                #Path line
                __, sv_region, nodes, nodes_len, __ = line.split("\t")

                if sv_region.startswith("_alt_"):
                    #Alternative path
                    continue
                
                nodes = extract_nodes(nodes)
                # if len(sv_region.split("_")) > 3:
                dict_of_sv_nodes = fill_sv_nodes(dict_of_sv_nodes, sv_region, nodes, dict_of_nodes_len)


    print(str(len(dict_of_sv_nodes)), "SVs")


    #2. Reading alignments and filling dictionary of informative alignments

    print("Reading gaf file...")

    dict_of_informative_aln = dict()
    list_of_identities = list()
    # dict = {sv_id : [[aln on ref allele], [aln on alt allele]]}

    ################## Tests
    # aln_reads_before_filters = dict()
    # aln_reads_after_filters = dict()
    ##################

    with open(gaf_file) as aln_file:
        for aln in aln_file:

            sv_based_alns = read_gaf_aln(aln, dict_of_sv_nodes, dict_of_nodes_len) #list of sv based aln(s) from graph-based aln

            for sv_based_aln in sv_based_alns:
                read, r_len, r_start, r_end, strand, sv_id, allele, a_len, a_start, a_end, a_identity = sv_based_aln

                ################## Tests
                # if read not in aln_reads_before_filters:
                #     # read_id = read.split(";")[0].split("Read=")[1]
                #     aln_reads_before_filters[read] = []
                # aln_reads_before_filters[read].append("\t".join([read, str(r_len), str(r_start), str(r_end), strand, sv_id, str(allele), str(a_len), str(a_start), str(a_end)]))
                ##################

                if allele is None:
                    continue

                #Breakpoint overlap filter
                if do_overlap_breakpoint(a_len, a_start, a_end, l_adj, d_over):
                    
                    #Semi-globality filter

                    if is_semiglobal_aln(r_len, r_start, r_end, a_len, a_start, a_end, d_end, l_adj):

                    #Ambiguous read alignment filter (remove redundant reads aligning on complementary ref)

                    # ... 
                    # find replacement score for mapping quality

                        if sv_id not in dict_of_informative_aln.keys():
                            dict_of_informative_aln[sv_id] = [[], []]

                        dict_of_informative_aln[sv_id][allele].append(aln.split("cg:Z:")[0])
                        list_of_identities.append(a_identity)

                        ################## Tests
                        # if read not in aln_reads_after_filters:
                        #     # read_id = read.split(";")[0].split("Read=")[1]
                        #     aln_reads_after_filters[read] = []
                        # aln_reads_after_filters[read].append("\t".join([read, str(r_len), str(r_start), str(r_end), strand, sv_id, str(allele), str(a_len), str(a_start), str(a_end)]))
                        ##################
    

    #3. Curating informative alignments on identity %

    min_id_value = round(stats.quantiles(list_of_identities, n=4)[0], 3)

    for sv, allele_alns in dict_of_informative_aln.items():
        for allele in [0, 1]:

            curated_alns = []
            for aln in allele_alns[allele]:
                identity = float(aln.rstrip().split("\t")[14].split(":")[-1])

                if identity > min_id_value:
                    curated_alns.append(aln)

            dict_of_informative_aln[sv][allele] = curated_alns


    ### Stats
    aln_nb = 0
    for sv in dict_of_informative_aln:
        aln_nb += len(dict_of_informative_aln[sv][0]) + len(dict_of_informative_aln[sv][1])
    
    #Save dict
    with open(output_aln_dict, 'w') as file:
        file.write(json.dumps(dict_of_informative_aln, sort_keys=True, indent=4))
    print(str(aln_nb), "informative aln saved")

    ################## Tests
    # aln_reads_dict = {"before_filters" : aln_reads_before_filters, "after_filters" : aln_reads_after_filters}
    # with open('reads_aln.json', 'w') as read_stats_file:
    #     read_stats_file.write(json.dumps(aln_reads_after_filters, sort_keys=True, indent=4))
    ##################  



def is_semiglobal_aln(r_len, r_start, r_end, a_len, a_start, a_end, d_end, l_adj):

    return any([(a_start <= d_end) and (a_len - d_end <= a_end), #allele inside read
                (r_start <= d_end) and (r_len - d_end <= r_end), #read inside allele
                (r_start <= d_end) and (a_len - d_end <= a_end), #read left aligned
                (a_start <= d_end) and (r_len - d_end <= r_end)]) #read right aligned
     
def do_overlap_breakpoint(a_len, start, end, l_adj, d_over):

    if (start + d_over) < l_adj < (end - d_over) or (start + d_over) < (a_len - l_adj) < (end - d_over):
        return True
    
    return False


def read_gaf_aln(gaf_line, sv_dict, nodes_len_dict):
    #Traduit l'alignement => path devient allele_svid ; p_len devient allele_len ; p_start devient allele_start ; p_end devient allele_end
    # => trouve le sv puis l'allèle puis convertis les positions

    read, r_len, r_start, r_end, strand, path, __, p_start, p_end, *__ = gaf_line.rstrip().split("\t")
    identity = float(gaf_line.rstrip().split("\t")[14].split(":")[-1])

    #Get nodes of path
    aln_nodes = extract_nodes(path, ordered=True) ## , strand
    if aln_nodes is None: #pb in path orientation
        return []
    elif len(aln_nodes) == 1: #aln on only 1 node can't overlap breakpt
        return []

    #Get list of sv the read mapped on
    target_svs = find_sv(aln_nodes, sv_dict)
    if len(target_svs) == 0:
        print('Did not find SV for aln:', gaf_line, 'with path:', path)
        return []

    #Translate graph-based aln to sv-based aln
    sv_based_alns = list()
    for sv_id in target_svs:
        sv_type = sv_id.split("_")[0]

        #DELETIONS
        if sv_type == 'DEL':
            allele = find_DEL_allele(aln_nodes, sv_id, sv_dict)
            a_len, a_start, a_end = get_aln_pos_on_allele_DEL(sv_id, allele, aln_nodes, int(p_start), int(p_end), sv_dict, nodes_len_dict)

        sv_id = "_".join(sv_id.split("_")[1:]) #remove sv type in sv id
        sv_based_alns.append([read, int(r_len), int(r_start), int(r_end), strand, sv_id, allele, a_len, a_start, a_end, identity])

    return sv_based_alns

def find_sv(aln_nodes, sv_dict):
    """
    Find SV from node in dict of SV nodes
    """
    target_sv = set()
    for sv in sv_dict:
        if aln_nodes[1] in sv_dict[sv][0]:
            target_sv.add(sv)
    return list(target_sv)

def find_DEL_allele(nodes, sv, sv_dict):
    #Look for ref specific nodes in aln nodes
    for ref_specific_node in sv_dict[sv][1]: #sv nodes = ref specific nodes
        if ref_specific_node in nodes:
            return 0

    #Look for alt specific node succession in aln nodes
    left_breakpt_node = str(int(sv_dict[sv][1][0]) - 1) #node preceeding first ref specific node
    right_breakpt_node = str(int(sv_dict[sv][1][-1]) + 1) #node following last ref specific node
    alt_path_node_succession = "-".join([left_breakpt_node, right_breakpt_node])
    aln_node_succession = "-".join(list(map(str, nodes)))

    if alt_path_node_succession in aln_node_succession:
        return 1
    
    #Aln do not inform on sv allele
    return None

def get_aln_pos_on_allele_DEL(sv, allele, aln_nodes, start, end, sv_dict, node_len_dict):

    if allele is None:
        return None, None, None

    a_pos_on_first_node = sv_dict[sv][2][0] #allele start on sv region's first node
    a_end_on_last_node = sv_dict[sv][2][1] #allele end on sv region's last node
    a_len = 0 #allele length
    a_start = None #aln start on allele
    first_aln_node_on_allele = 0 #first aln node on allele

    while (a_start is None) and (first_aln_node_on_allele < len(aln_nodes)):

        for sv_node in sv_dict[sv][0]:
            #Look for aln start node
            if sv_node == aln_nodes[first_aln_node_on_allele]:
                a_start = a_len - a_pos_on_first_node + start - 1 #aln start on allele
                a_end = a_len - a_pos_on_first_node + end - 1 #aln end on allele

            #Add nodes length to allele length
            if (allele == 1) and (sv_node in sv_dict[sv][1]):
                #Don't count ref nodes length in alt allele length
                continue

            else:
                if sv_node == sv_dict[sv][0][-1]:
                    #Add only remaining allele length for last node of sv region
                    a_len += a_end_on_last_node + 1
                else:
                    a_len += node_len_dict[sv_node]
        
        first_aln_node_on_allele += 1

    return a_len, a_start, a_end

def fill_sv_nodes(node_dict, sv_region, nodes, dict_of_nodes_len):
    l_adj = 5000
    chrom = sv_region.split("_")[1]
    first_sv_start = int(sv_region.split('_')[2].split('-')[1]) #start of first sv on chrom
    sv_number = len(sv_region.split("_")) - 2

    for sv in range(1, sv_number + 1):
        sv_type, pos, length = sv_region.split('_')[1 + sv].split('-')
        sv_id = "_".join([sv_type, chrom, "-".join([pos, length])])
        node_dict[sv_id] = [[], [], []] #[[sv region nodes], [sv nodes], [sv region start on first region node, sv region end on last region node]]

        #Get coordinates of sv and sv region
        if sv_type == 'DEL':
            region_start = int(pos) - first_sv_start #start of current sv REGION on graph
            sv_start = int(pos) - first_sv_start + l_adj #start of current sv on graph
            if int(length) <= 2 * l_adj:
                sv_end = sv_start + int(length) - 1 #end of current sv on graph
            else:
                sv_end = sv_start + 2 * l_adj - 1 #end if sv sequence cut
            region_end = sv_end + l_adj #end of current sv REGION on graph
        
        # elif sv_type == 'ins':

        #Get nodes of sv
        pos_on_graph = 0
        for node in nodes:

            #Find first node of sv
            if (pos_on_graph < region_start + 1) and (pos_on_graph + dict_of_nodes_len[node] > region_start): #node contains sv region start
                node_dict[sv_id][0].append(node) #add node to sv region nodes

                region_start_on_node = region_start - pos_on_graph
                node_dict[sv_id][2].append(region_start_on_node)

                pos_on_graph += dict_of_nodes_len[node]

                while (pos_on_graph < region_end) and (node < nodes[-1]): #sv region end after start of node
                    node += 1
                    node_dict[sv_id][0].append(node)

                    #Get sv specific nodes
                    if (pos_on_graph > sv_start - 1) and (pos_on_graph < sv_end + 1): #node contains sv
                        node_dict[sv_id][1].append(node) #add node to sv nodes

                    pos_on_graph += dict_of_nodes_len[node]
                
                #Found all sv region nodes
                region_end_on_node = region_end - (pos_on_graph - dict_of_nodes_len[node - 1]) #DEL, INS, INV
                node_dict[sv_id][2].append(region_end_on_node)
                break

            pos_on_graph += dict_of_nodes_len[node]
    
    return node_dict

def fill_allele_specific_nodes(alt_path_line):
    pass #implement as function later #or not

def check_path_orientation(path, nodes):
    nodes = list(map(str, nodes))

    for i in range(len(nodes)):

        if i == 0:
            orient = path.split(nodes[i])[0]

        else:
            next_orient = path.split(nodes[i-1])[1].split(nodes[i])[0]
            if next_orient != orient:
                #Ambiguous path orientation
                # print("Warning: ambiguous orientation in alignment path", path, "- will be ignored")
                return ""

    return orient

def extract_nodes(str_node_list, *args, **kwargs):
    """
    Return nodes in list format from nodes in str format. 
    Can take path orientation into account if input is gaf path.

    Args:
        str_node_list (str) : node list in str, with nodes separated by special characters
        ordered (bool) : take path orientation into account or not

    Ex: '1+,2+,3+,4+' => [1, 2, 3, 4]   (path in gfa)
           '>1>2>3>4' => [1, 2, 3, 4]   (path in gaf)
    """
    check_orient = kwargs.get('ordered', False)

    nodes = re.split(r'\W+', str_node_list)

    if nodes[0] == "":
        nodes = nodes[1:]
    if nodes[-1] == "":
        nodes = nodes[:-1]
    
    nodes = list(map(int, nodes))

    if check_orient:
        orient = check_path_orientation(str_node_list, nodes)

        if (orient == "<") and (len(nodes) > 1):
            # strand = '-'
            nodes.reverse()
            return nodes #can return strand if needed
        
        elif (orient == ">"):
            # strand = '+'
            return nodes #can return strand if needed

        elif orient == "":
            #Multiple aln orientations in path
            print("Multiple path orientations not supported yet (function extract_nodes())")
            return None

    return nodes


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])