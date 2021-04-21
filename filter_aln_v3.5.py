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
        "-n", 
        "--sv_nodes", metavar="<sv_node_dict>", nargs=1, 
        help="sv nodes info", 
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
        "-v", 
        "--version", 
        metavar="<versionID", 
        type=str,
        required=False)

    args = parser.parse_args()
    
    if args.sv_nodes is not None:
        pass #load dict_of_DEL_sv_nodes from json

    if args.version:
        output_aln_dict = 'informative_aln_v3.5_' + args.version + '.json'
    else:
        output_aln_dict = 'informative_aln_v3.5.json'

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
    #dict = {sv_id : [[sv region nodes], [sv nodes], [sv region start on first region node, sv region end on last region node], [first node of component, last node of component]]}
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
                read, r_len, r_start, r_end, strand, sv_id, allele, a_len, a_start, a_end, a_coords_for_semiglob, a_identity = sv_based_aln

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

                    if is_semiglobal_aln(r_len, r_start, r_end, a_coords_for_semiglob, dict_of_sv_nodes["DEL_" + sv_id][3], d_end, dict_of_nodes_len):

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



def is_semiglobal_aln(r_len, r_start, r_end, a_coords_for_semiglob, list_component_first_last_node, d_end, nodes_len_dict):

    # a_coords_for_semiglob = [[a_start_node, a_start_on_node], [a_end_node, a_end_on_node]]
    # list_component_first_last_node = [component first node, component last node]

    global_on_component_left = (list_component_first_last_node[0] == a_coords_for_semiglob[0][0]) and (a_coords_for_semiglob[0][1] <= d_end)
    global_on_component_right = (list_component_first_last_node[1] == a_coords_for_semiglob[1][0]) and (nodes_len_dict[a_coords_for_semiglob[1][0]] - d_end <= a_coords_for_semiglob[1][1])

    return any([(global_on_component_left) and (global_on_component_right), #allele inside read
                (r_start <= d_end) and (r_len - d_end <= r_end), #read inside allele
                (r_start <= d_end) and (global_on_component_right), #read left aligned
                (global_on_component_left) and (r_len - d_end <= r_end)]) #read right aligned
     
def do_overlap_breakpoint(a_len, start, end, l_adj, d_over):

    if (start + d_over) < l_adj < (end - d_over) or (start + d_over) < (a_len - l_adj) < (end - d_over):
        return True
    
    return False


def read_gaf_aln(gaf_line, sv_dict, nodes_len_dict):
    #Traduit l'alignement => path devient allele_svid ; p_len devient allele_len ; p_start devient allele_start ; p_end devient allele_end
    # => trouve le sv puis l'allèle puis convertis les positions

    read, r_len, r_start, r_end, strand, path, p_len, p_start, p_end, *__ = gaf_line.rstrip().split("\t")
    identity = float(gaf_line.rstrip().split("\t")[14].split(":")[-1])

    #Get nodes of path
    aln_nodes = extract_nodes(path, ordered=True) ## , strand
    if aln_nodes is None: #pb in path orientation
        return []
    elif len(aln_nodes) == 1: #aln on only 1 node can't pass breakpt overap filter
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

        #DELETIONS, INSERTIONS
        if sv_type in ["DEL", "INS"]:
            allele = find_DEL_INS_allele(aln_nodes, sv_id, sv_dict, sv_type)
            a_len, a_start, a_end, a_coords_for_semiglob = get_aln_pos_on_DEL_INS_allele(sv_id, sv_type, allele, aln_nodes, int(p_len), int(p_start), int(p_end), sv_dict, nodes_len_dict)

        #INVERSIONS
        if sv_type == "INV":
            continue
            #ALLELE : 1.determine global path orientation, 2.if global orient = "f": if ">first sv node" or ">last sv node": allele 0; 
            #POS ON ALLELE : same as DEL/INS but if aln starts/ends on sv node and allele 1

        sv_id = "_".join(sv_id.split("_")[1:]) #remove sv type in sv id
        sv_based_alns.append([read, int(r_len), int(r_start), int(r_end), strand, sv_id, allele, a_len, a_start, a_end, a_coords_for_semiglob, identity])

    return sv_based_alns

def find_sv(aln_nodes, sv_dict):
    """
    Find SV from node in dict of SV nodes
    """
    target_sv = set()
    for sv in sv_dict:
        for i in range(1, len(aln_nodes), 2): #only search matching sv for every second node of aln path
            if aln_nodes[i] in sv_dict[sv][0]:
                target_sv.add(sv)
    return list(target_sv)

def find_DEL_INS_allele(aln_nodes, sv, sv_dict, sv_type):
    '''Works for DEL and INS types'''
    #Look for sv specific nodes in aln nodes
    for sv_node in sv_dict[sv][1]: #sv nodes = sv sequence nodes
        if sv_node in aln_nodes:
            if sv_type == "DEL":
                return 0
            elif sv_type == "INS":
                return 1

    #Look for sv flanking nodes succession in aln nodes
    left_breakpt_node = str(int(sv_dict[sv][1][0]) - 1) #node preceeding first sv specific node
    right_breakpt_node = str(int(sv_dict[sv][1][-1]) + 1) #node following last sv specific node
    flanking_nodes_succession = "-".join([left_breakpt_node, right_breakpt_node])
    aln_node_succession = "-".join(list(map(str, aln_nodes)))

    if flanking_nodes_succession in aln_node_succession:
        if sv_type == "DEL":
            return 1
        elif sv_type == "INS":
            return 0
    
    #Aln do not inform on sv allele
    return None

def get_aln_pos_on_DEL_INS_allele(sv, sv_type, allele, aln_nodes, length, start, end, sv_dict, node_len_dict):

    if allele is None:
        return None, None, None, [None]

    ##########
    a_pos_on_first_node = sv_dict[sv][2][0] #allele start on sv region's first node
    a_end_on_last_node = sv_dict[sv][2][1] #allele end on sv region's last node
    a_start = None #aln start on allele
    first_aln_node_on_allele = 0 #first aln node on allele

    while (a_start is None) and (first_aln_node_on_allele < len(aln_nodes)):
        a_len = 0 #allele length

        for sv_node in sv_dict[sv][0]:
            #Look for aln start node
            if sv_node == aln_nodes[first_aln_node_on_allele]:
                a_start = a_len - a_pos_on_first_node + start - 1 #aln start on allele
                a_end = a_len - a_pos_on_first_node + end - 1 #aln end on allele

            #Add nodes length to allele length
            if (sv_type == "DEL") and (allele == 1) and (sv_node in sv_dict[sv][1]):
                #Don't count ref nodes length in alt allele length
                continue
            elif (sv_type == "INS") and (allele == 0) and (sv_node in sv_dict[sv][1]):
                #Don't count alt nodes length in ref allele length
                continue

            else:
                if sv_node == sv_dict[sv][0][-1]:
                    #Add only remaining allele length for last node of sv region
                    a_len += a_end_on_last_node + 1
                else:
                    a_len += node_len_dict[sv_node]
        
        first_aln_node_on_allele += 1
        ##########

    if sv_type in ["DEL", "INS"]:
        a_start_node = aln_nodes[0]
        a_start_on_node = start

        a_end_node = aln_nodes[-1]
        a_end_on_node = length - end

    return a_len, a_start, a_end, [[a_start_node, a_start_on_node], [a_end_node, a_end_on_node]]

def fill_sv_nodes(node_dict, sv_region, nodes, dict_of_nodes_len):
    l_adj = 5000
    chrom = sv_region.split("_")[1]
    first_sv_start = int(sv_region.split('_')[2].split('-')[1]) #start of first sv on chrom
    sv_number = len(sv_region.split("_")) - 2

    for sv in range(1, sv_number + 1):
        sv_type, pos, length = sv_region.split('_')[1 + sv].split('-')
        sv_id = "_".join([sv_type, chrom, "-".join([pos, length])])
        node_dict[sv_id] = [[], [], [], []] #[[sv region nodes], [sv nodes], [sv region start on first region node, sv region end on last region node], [component first node, component last node]]

        #Get coordinates of sv and sv region
        if sv_type in ["DEL", "INS", "INV"]:
            region_start = int(pos) - first_sv_start #start of current sv REGION on graph
            sv_start = int(pos) - first_sv_start + l_adj #start of current sv on graph
            if int(length) <= 2 * l_adj:
                sv_end = sv_start + int(length) - 1 #end of current sv on graph
            else:
                sv_end = sv_start + 2 * l_adj - 1 #end if sv sequence cut
            region_end = sv_end + l_adj #end of current sv REGION on graph

        #Get nodes of sv (works for DEL, INS, INV)
        pos_on_graph = 0
        for node in nodes:

            #Find first node of sv REGION
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
                
                #Found all sv region nodes (exited while loop)
                region_end_on_node = region_end - (pos_on_graph - dict_of_nodes_len[node - 1]) #works for DEL, INS, INV
                node_dict[sv_id][2].append(region_end_on_node)
                break

            pos_on_graph += dict_of_nodes_len[node]
        
        node_dict[sv_id][3] = [nodes[0], nodes[-1]]
    
    return node_dict

def fill_allele_specific_nodes(alt_path_line):
    pass #implement as function later #or not

def check_global_path_orientation(gaf_path, nodes):

    if (nodes[0] > nodes[-1]) and (gaf_path[0] == "<"):
        orient = "r" #reverse
    else:
        orient = "f" #forward

    return orient

def extract_nodes(path, *args, **kwargs):
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

    nodes = re.split(r'\W+', path)

    if nodes[0] == "":
        nodes = nodes[1:]
    if nodes[-1] == "":
        nodes = nodes[:-1]
    
    nodes = list(map(int, nodes))

    if check_orient:
        orient = check_global_path_orientation(path, nodes)

        if (orient == "r") and (len(nodes) > 1):
            # strand = '-'
            nodes.reverse()
            return nodes #can return strand if needed
        
        elif (orient == "f"):
            # strand = '+'
            return nodes #can return strand if needed

    return nodes


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])