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
        "-L", "--ladj", metavar="<len_adjacent_seq>", nargs=1, required=False, type=int, default=5000
    )
    parser.add_argument(
        "-O", "--dover", metavar="<min_breakpoint_overlap>", nargs=1, required=False, type=int, default=[100]
    )
    parser.add_argument(
        "-E", "--dend", metavar="<max_dist_semiglobal>", nargs=1, required=False, type=int, default=[100]
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
    d_over = args.dover[0]
    d_end = args.dend[0]


    #1. Preparing dict of SV graph info
    print("Reading gfa file...")

    if args.sv_nodes is not None:
        pass #don't recreate dict (real time save or not?)

    dict_of_sv_graph_info = dict()
    #dict = {sv_id : [[ref breakpt edges], [alt breakpt edges], [first node of graph, last node of graph]]}
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
                __, sv_region, nodes, __, __ = line.split("\t")

                if sv_region.startswith("_alt_"):
                    #Alternative path
                    continue
                
                nodes = extract_nodes(nodes)
                # if len(sv_region.split("_")) > 3:
                dict_of_sv_graph_info = add_region_svs(dict_of_sv_graph_info, sv_region, nodes, dict_of_nodes_len, l_adj)

                # else:
                #     sv_id = "_".join([sv_region.split("_")[1], sv_region.split("_")[2].split("-")[1:]])
                #     dict_of_sv_nodes[sv_id] = nodes

            # if line[0] == "L":
            #     #Link line
            #     continue

    #Save dict
    with open('sv_graph_info.json', 'w') as sv_dict_file:
        sv_dict_file.write(json.dumps(dict_of_sv_graph_info, sort_keys=True, indent=4))
    
    print(str(len(dict_of_sv_graph_info)), "SVs")


    #2. Reading alignments and filling dictionary of informative alignments

    print("Reading gaf file...")

    dict_of_informative_aln = dict()
    # dict = {sv_id : [[aln on ref allele], [aln on alt allele]]}
    list_of_identities = list()

    ################## Tests
    total_number_alns = 0
    alns_passing_breakpt = 0
    alns_passing_semiglobal = 0
    alns_passing_identity = 0
    ##################

    ################## Tests
    # aln_reads_before_filters = dict()
    # aln_reads_after_filters = dict()
    ##################

    with open(gaf_file) as aln_file:
        for aln in aln_file:

            # print("\n", aln, "\n")

            ################## Tests
            total_number_alns += 1
            ##################

            sv_based_alns = read_gaf_aln(aln, dict_of_sv_graph_info, dict_of_nodes_len) #list of sv based aln(s) from graph-based aln

            for sv_based_aln in sv_based_alns:
                read, read_len, start_on_read, end_on_read, sv_id, allele, aln_path, start_on_first_node, end_on_last_node, a_identity = sv_based_aln

                # print(sv_based_aln)

                ################## Tests
                # if read not in aln_reads_before_filters:
                #     # read_id = read.split(";")[0].split("Read=")[1]
                #     aln_reads_before_filters[read] = []
                # aln_reads_before_filters[read].append("\t".join([read, str(r_len), str(r_start), str(r_end), strand, sv_id, str(allele), str(a_len), str(a_start), str(a_end)]))
                ##################

                if allele is None:
                    continue

                #Breakpoint overlap filter

                #v2
                pass_breakpt_filter = False
                for breakpt_edge in dict_of_sv_graph_info[sv_id][allele]:
                    if do_overlap_breakpoint(breakpt_edge, aln_path, start_on_first_node, end_on_last_node, dict_of_nodes_len, d_over):
                        pass_breakpt_filter = True
                        break

                if pass_breakpt_filter:

                    ################## Tests
                    # print("Passed breakpt filter")
                    alns_passing_breakpt += 1
                    ##################
                
                #v1
                # if any([do_overlap_breakpoint(breakpt_edge, aln_path, start_on_first_node, end_on_last_node, dict_of_nodes_len, d_over)] for breakpt_edge in dict_of_sv_graph_info[sv_id][allele]):

                    #Semi-globality filter

                    if is_semiglobal_aln(read_len, start_on_read, end_on_read, dict_of_sv_graph_info[sv_id][2], aln_path, start_on_first_node, end_on_last_node, dict_of_nodes_len, d_end):

                        ################## Tests
                        # print("Passed semiglobal filter")
                        alns_passing_semiglobal += 1
                        ##################

                    #Ambiguous read alignment filter (remove redundant reads aligning on complementary ref)

                    # ... 
                    # find replacement score for mapping quality

                        # Modify sv_id to match with svjedi_genotype_prediction.py
                        # chrom_DEL_start-end -> chrom_start-length
                        sv_type, sv_chrom, sv_pos = sv_id.split("_")
                        sv_start, sv_end = sv_pos.split("-")
                        if sv_type == "DEL":
                            sv_length = str(int(sv_end) - int(sv_start))
                        sv_id = "_".join([sv_chrom, "-".join([sv_start, sv_length])])


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

            ################## Tests
            alns_passing_identity += len(curated_alns)
            ##################


    ### Stats
    aln_nb = 0
    for sv in dict_of_informative_aln:
        aln_nb += len(dict_of_informative_aln[sv][0]) + len(dict_of_informative_aln[sv][1])

    print(total_number_alns, "alns in gaf file.")
    print(alns_passing_breakpt, "alns passing breakpoint filter.")
    print(alns_passing_semiglobal, "alns passing semi-globality filter.")
    print(alns_passing_identity, "alns passing identity filter.")
    
    #Save dict
    with open(output_aln_dict, 'w') as file:
        file.write(json.dumps(dict_of_informative_aln, sort_keys=True, indent=4))
    print(str(aln_nb), "informative aln saved")

    ################## Tests
    # aln_reads_dict = {"before_filters" : aln_reads_before_filters, "after_filters" : aln_reads_after_filters}
    # with open('reads_aln.json', 'w') as read_stats_file:
    #     read_stats_file.write(json.dumps(aln_reads_after_filters, sort_keys=True, indent=4))
    ##################  



def is_semiglobal_aln(r_len, r_start, r_end, first_last_graph_nodes, aln_path, start_on_first_node, end_on_last_node, nodes_len_dict, d_end):
    '''
    Check if aln starts at the start (+/- d_end) of the read OR of the graph, AND ends at the end (+/- d_end) of the read OR of the graph.
    '''
    aln_first_node = int(aln_path.split(",")[0][0:-1])
    aln_last_node = int(aln_path.split(",")[-1][0:-1])

    # global_on_read_start = (r_start <= d_end)
    # global_on_read_end = (r_len - d_end <= r_end)

    # global_on_graph_start = (first_aln_node == first_last_graph_nodes[0] and start_on_first_node <= d_end)
    # global_on_graph_end = (last_aln_node == first_last_graph_nodes[1] and nodes_len_dict[last_aln_node] - d_end <= end_on_last_node)

    # if global_on_graph_start or global_on_graph_end:
    #     print(global_on_read_start, "|", global_on_graph_start, "&", global_on_read_end, "|", global_on_graph_end)

    #v4
    semiglobal_on_start = (r_start <= d_end) or (aln_first_node == first_last_graph_nodes[0] and start_on_first_node <= d_end)
    semiglobal_on_end = (r_len - d_end <= r_end) or (aln_last_node == first_last_graph_nodes[1] and nodes_len_dict[aln_last_node] - d_end <= end_on_last_node)
    #v5
    # semiglobal_on_start = (r_start <= d_end) or (start_on_first_node <= d_end)
    # semiglobal_on_end = (r_len - d_end <= r_end) or (nodes_len_dict[aln_last_node] - d_end <= end_on_last_node)

    return (semiglobal_on_start and semiglobal_on_end)
     
def do_overlap_breakpoint(breakpt_edge, aln_path, start_on_first_node, end_on_last_node, nodes_len_dict, d_over):

    if breakpt_edge in aln_path:

        breakpt_left_node = int(breakpt_edge.split(",")[0][0:-1])
        breakpt_right_node = int(breakpt_edge.split(",")[1][0:-1])
        aln_first_node = int(aln_path.split(",")[0][0:-1])
        aln_last_node = int(aln_path.split(",")[-1][0:-1])

        overlap_on_left = (breakpt_left_node > aln_first_node) or (nodes_len_dict[aln_first_node] - start_on_first_node >= d_over) #breakpt_left_node either > or = aln_first_node
        overlap_on_right = (breakpt_right_node < aln_last_node) or (end_on_last_node + 1 >= d_over) #breakpt_right_node either < or = aln_last_node

        return overlap_on_left and overlap_on_right
    
    return False


def read_gaf_aln(gaf_line, sv_graphInfo_dict, nodes_len_dict):
    #Traduit l'alignement => path devient allele_svid ; p_len devient allele_len ; p_start devient allele_start ; p_end devient allele_end
    # => trouve le sv puis l'allèle puis convertis les positions

    read, read_len, start_on_read, end_on_read, __, path, __, p_start, p_end, *__ = gaf_line.rstrip().split("\t")
    identity = float(gaf_line.rstrip().split("\t")[14].split(":")[-1])

    #Translate aln path
    aln_nodes = extract_nodes(path)
    aln_path = translate_path(path, aln_nodes)
    if len(aln_path.split(",")) == 1: #aln on only 1 node can't overlap breakpt on both sides
        return []

    #Get list of sv the read mapped on
    target_sv_allele = find_sv_allele(aln_path, sv_graphInfo_dict)
    if len(target_sv_allele) == 0:
        #print('Did not find SV for aln:', gaf_line, 'with path:', path)
        return []

    #Translate graph-based aln to sv-based aln
    sv_based_alns = list()
    for (sv_id, allele) in target_sv_allele:
        
        #sv_type = sv_id.split("_")[0]
        #DELETIONS
        # if sv_type == 'DEL':
        #     a_len, a_start, a_end = get_aln_pos_on_allele_DEL(sv_id, allele, aln_path, int(p_start), int(p_end), sv_breakpt_dict, nodes_len_dict)

        start_on_first_node = int(p_start)
        end_on_last_node = int(p_end)
        for node in aln_nodes[:-1]:
            end_on_last_node -= nodes_len_dict[node]

        #sv_id = "_".join(sv_id.split("_")[1:]) #remove sv type in sv id
        sv_based_alns.append([read, int(read_len), int(start_on_read), int(end_on_read), sv_id, allele, aln_path, start_on_first_node, end_on_last_node, identity])

    return sv_based_alns

def find_sv_allele(aln_path, sv_graphInfo_dict):
    """
    Find SV from path in dict of SV breakpt edges
    """
    target_sv_allele = set()
    for sv in sv_graphInfo_dict.keys(): #{sv : [[ref breakpt edges], [alt breakpt edges], [first and last nodes of region graph]]}

        for ref_breakpt_edge in sv_graphInfo_dict[sv][0]:
            if ref_breakpt_edge in aln_path:
                target_sv_allele.add((sv, 0))

        for alt_breakpt_edge in sv_graphInfo_dict[sv][1]:
            if alt_breakpt_edge in aln_path:
                target_sv_allele.add((sv, 1))
        
        if ((sv, 0) in target_sv_allele) and ((sv, 1) in target_sv_allele): #ambiguous aln on sv
            target_sv_allele.remove((sv, 0))
            target_sv_allele.remove((sv, 1))

    return list(target_sv_allele)

def add_region_svs(sv_graphInfo_dict, sv_region, region_nodes, dict_of_nodes_len, l_adj):
    chrom = sv_region.split("_")[1]
    first_sv_start = int(sv_region.split('_')[2].split('-')[1]) #start of first sv on chrom
    sv_number = len(sv_region.split("_")) - 2

    #Get info for each SV in the graph region
    for sv in range(1, sv_number + 1):
        sv_type, pos, end = sv_region.split('_')[1 + sv].split('-')
        sv_id = "_".join([sv_type, chrom, "-".join([pos, end])])
        sv_graphInfo_dict[sv_id] = [ [], [], [] ] #[ [ref breakpt edges], [alt breakpt edges],  [first node of graph, last node of graph]] for DEL, INS, INV

        #Get coordinates of sv and sv region
        sv_start = int(pos) - first_sv_start + l_adj #start of current sv on graph
        sv_end = int(end) - first_sv_start + l_adj #end of current sv on graph

        # if sv_type == 'DEL' or sv_type == 'INS':
        #     if int(length) <= 2 * l_adj:
        #         sv_end = sv_start + int(length) - 1 #end of current sv on graph
        #     else:
        #         sv_end = sv_start + 2 * l_adj - 1 #end if sv sequence cut
        
        # elif sv_type == 'ins':

        #Get breakpoint nodes of sv (works for DEL, INS, INV)
        pos_on_graph = 0
        for node in region_nodes:

            if pos_on_graph == sv_start: #start of node == start of sv
                longAllele_left_breakpt_edge = ",".join([str(node - 1) + "+",  str(node) + "+"]) #ex : '1+,2+' (node before breakpt, node after breakpt)
            
            if (pos_on_graph + dict_of_nodes_len[node]) == sv_end: #end of node == end of sv
                longAllele_right_breakpt_edge = ",".join([str(node) + "+",  str(node + 1) + "+"])

            pos_on_graph += dict_of_nodes_len[node]
        
        #Fill dict with sv breakpt edges and first and last nodes of region graph
        if sv_type == 'DEL':
            ref_left_breakpt_edge = longAllele_left_breakpt_edge
            ref_right_breakpt_edge = longAllele_right_breakpt_edge
            alt_breakpt_edge = ",".join([longAllele_left_breakpt_edge.split(",")[0], longAllele_right_breakpt_edge.split(",")[1]])

            sv_graphInfo_dict[sv_id] = [[ref_left_breakpt_edge, ref_right_breakpt_edge], [alt_breakpt_edge], [region_nodes[0], region_nodes[-1]]]
        
        elif sv_type == 'INS':
            ref_breakpt_edge = ",".join([longAllele_left_breakpt_edge.split(",")[0], longAllele_right_breakpt_edge.split(",")[1]])
            alt_left_breakpt_edge = longAllele_left_breakpt_edge
            alt_right_breakpt_edge = longAllele_right_breakpt_edge

            sv_graphInfo_dict[sv_id] = [[ref_breakpt_edge], [alt_left_breakpt_edge, alt_right_breakpt_edge], [region_nodes[0], region_nodes[-1]]]
        
        elif sv_type == 'INV':
            ref_left_breakpt_edge = longAllele_left_breakpt_edge
            ref_right_breakpt_edge = longAllele_right_breakpt_edge
            inverted_left_breakpt_node = longAllele_left_breakpt_edge.split(",")[1].replace("+", "-") #left breakpt = '1+,2+' => inverted node = '2-'
            inverted_right_breakpt_node = longAllele_right_breakpt_edge.split(",")[0].replace("+", "-") #right breakpt = '3+,4+' => inverted node = '3-'
            alt_left_breakpt_edge = ",".join([longAllele_left_breakpt_edge.split(",")[0], inverted_left_breakpt_node])
            alt_right_breakpt_edge = ",".join([inverted_right_breakpt_node, longAllele_right_breakpt_edge.split(",")[1]])

            sv_graphInfo_dict[sv_id] = [[ref_left_breakpt_edge, ref_right_breakpt_edge], [alt_left_breakpt_edge, alt_right_breakpt_edge], [region_nodes[0], region_nodes[-1]]]
    
    return sv_graphInfo_dict

def extract_nodes(str_node_list):
    """
    Return nodes in list format from nodes in str format. 

    Args:
        str_node_list (str) : node list in str, with nodes separated by special characters

    Ex: '1+,2+,3+,4+' => [1, 2, 3, 4]   (path in gfa)
           '>1>2>3>4' => [1, 2, 3, 4]   (path in gaf)
    """

    nodes = re.split(r'\W+', str_node_list)

    if nodes[0] == "":
        nodes = nodes[1:]
    if nodes[-1] == "":
        nodes = nodes[:-1]
    
    nodes = list(map(int, nodes))
    return nodes

def translate_path(gaf_path, extracted_nodes):

    strand = {"for" : {">" : "+", "<" : "-"}, "rev" : {">" : "-", "<" : "+"}}
    
    translated_path = ''

    #Check for orientation of aln
    if (extracted_nodes[0] > extracted_nodes[-1]) and (gaf_path[0] == "<"):
        aln_orient = "rev"
        extracted_nodes.reverse()
    else:
        aln_orient = "for"

    for node in extracted_nodes:
        path_orientation_on_node = gaf_path.split(str(node))[0][-1]

        if translated_path == '':
            translated_path = str(node) + strand[aln_orient][path_orientation_on_node]
        else:
            translated_path = ",".join([translated_path, str(node) + strand[aln_orient][path_orientation_on_node]])
    
    return translated_path


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])