#!/usr/bin/env python3

import sys
import argparse
import re
import json

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

    # parser.add_argument(
    #     "-v", "--vcf", metavar="<vcf_file>", nargs=1, help="vcf file", required=True
    # )

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

    args = parser.parse_args()
    
    if args.sv_nodes is not None:
        pass #load dict_of_DEL_sv_nodes from json

    # vcf_file = args.vcf[0]
    gfa_file = args.gfa[0]
    gaf_file = args.gaf[0]
    l_adj = args.ladj
    d_over = args.dover
    d_end = args.dend

    #1. Preparing dict of SV nodes
    print("Reading gfa file...")

    if args.sv_nodes is not None:
        pass #don't recreate dict (real time save or not?)

    dict_of_DEL_sv_nodes = dict()
    #dict = {sv_id : [[nodes], [ref specific node], [alt specific node]]}
    dict_of_nodes_len = dict()
    #dict = {node_id : node_len}
    
    #1.1. Writing node list for each SV 

    with open(gfa_file) as graph_file:
        for line in graph_file:  

            if line[0] == "S":
                #Node line
                __, node_id, node_seq = line.rstrip().split("\t")
                dict_of_nodes_len[node_id] = len(node_seq)

            if line[0] == "P":
                #Path line
                __, sv_region, nodes, nodes_len, __ = line.split("\t")

                if sv_region.startswith("_alt_"):
                    #Alternative path
                    continue
                
                ########## def fill_sv_nodes()
                # __, sv_id, sv_type, sv_len = sv_region.split("_")
                sv_id = "_".join(sv_region.split("_")[1:])
                    
                # if sv_type == "del":
                    # dict_of_DEL_sv_nodes[sv_id] = [[], [], []]
                    # dict_of_DEL_sv_nodes[sv_id][0] = extract_nodes(nodes)

                dict_of_DEL_sv_nodes[sv_id] = [[], [], []]
                dict_of_DEL_sv_nodes[sv_id][0] = extract_nodes(nodes)
                ##########

            # if line[0] == "L":
            #     #Link line
            #     continue

    #1.2. Writing allele specific nodes for each SV

    with open(gfa_file) as graph_file:
        for line in graph_file:

            if line[0] == "P":
                #Path line
                __, sv_region, nodes, nodes_len, __ = line.split("\t")

                if sv_region.startswith("_alt_"):
                    #Alternative path

                    ########## def fill_allele_specific_nodes()
                    __, __, __, allele = sv_region.split("_")

                    if (allele == "0") and (nodes != ""):
                        #Ref specific nodes of DELETION

                        nodes = extract_nodes(nodes)
                        sv_id = find_sv(nodes[0], dict_of_DEL_sv_nodes)
                        dict_of_DEL_sv_nodes[sv_id][1] = nodes
                    
                    # elif (allele == "1") and (nodes != ""):
                    #     #Alt specific nodes of INSERTION
                    #     dict_of_INS_sv_nodes[sv_id][2] = extract_nodes(nodes)
                    ##########

    #Save dict
    dict_of_sv_nodes = dict()
    dict_of_sv_nodes['DEL'] = dict_of_DEL_sv_nodes
    # with open('sv_nodes_info.json', 'w') as sv_dict_file:
    #     sv_dict_file.write(json.dumps(dict_of_sv_nodes, sort_keys=True, indent=4))
    
    print(str(len(dict_of_sv_nodes['DEL'])), "deletion SVs")



    #2. Reading alignments

    print("Reading gaf file...")

    dict_of_informative_aln = dict()
    # dict = {sv_id : [[aln on ref allele], [aln on alt allele]]}

    ################## Tests
    # aln_reads_before_filters = dict()
    # aln_reads_after_filters = dict()
    ##################

    with open(gaf_file) as aln_file:
        for aln in aln_file:

            read, r_len, r_start, r_end, strand, sv_id, allele, a_len, a_start, a_end = read_gaf_aln(aln, dict_of_sv_nodes, dict_of_nodes_len)

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

                    ################## Tests
                    # if read not in aln_reads_after_filters:
                    #     # read_id = read.split(";")[0].split("Read=")[1]
                    #     aln_reads_after_filters[read] = []
                    # aln_reads_after_filters[read].append("\t".join([read, str(r_len), str(r_start), str(r_end), strand, sv_id, str(allele), str(a_len), str(a_start), str(a_end)]))
                    ##################
    
    ### Stats
    aln_nb = 0
    for sv in dict_of_informative_aln:
        aln_nb += len(dict_of_informative_aln[sv][0]) + len(dict_of_informative_aln[sv][1])
    
    #Save dict
    with open('informative_aln.json', 'w') as file:
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


def read_gaf_aln(gaf_line, all_sv_dict, nodes_len_dict):
    #Traduit l'alignement => path devient allele_svid ; p_len devient allele_len ; p_start devient allele_start ; p_end devient allele_end
    # => trouve le sv puis l'allèle puis convertis les positions

    read, r_len, r_start, r_end, strand, path, p_len, p_start, p_end, *__ = gaf_line.rstrip().split("\t")

    aln_nodes = extract_nodes(path, ordered=True) ## , strand

    #DELETIONS
    sv_id = find_sv(aln_nodes[0], all_sv_dict['DEL'])
    allele = find_DEL_allele(aln_nodes, sv_id, all_sv_dict['DEL'])
    a_len, a_start, a_end = get_aln_pos_on_allele_DEL(sv_id, allele, aln_nodes, int(p_start), int(p_end), all_sv_dict['DEL'], nodes_len_dict)

    return read, int(r_len), int(r_start), int(r_end), strand, sv_id, allele, a_len, a_start, a_end

def find_sv(node, sv_type_dict):
    """
    Find SV from node in dict of SV nodes
    """
    for sv in sv_type_dict:

        if node in sv_type_dict[sv][0]:
            return sv

def find_DEL_allele(nodes, sv, sv_del_dict):

    #Look for ref specific nodes in aln nodes
    for ref_specific_node in sv_del_dict[sv][1]:
        if ref_specific_node in nodes:
            return 0

    #Look for alt specific node succession in aln nodes
    left_breakpt_node = str(int(sv_del_dict[sv][1][0]) - 1) #node preceeding first ref specific node
    right_breakpt_node = str(int(sv_del_dict[sv][1][-1]) + 1) #node following last ref specific node
    alt_path_node_succession = "-".join([left_breakpt_node, right_breakpt_node])
    aln_nodes_succession = "-".join(nodes)

    if alt_path_node_succession in aln_nodes_succession:
        return 1
    
    #Read not aligned on allele specific node/path
    return None

def get_aln_pos_on_allele_DEL(sv, a, aln_nodes, start, end, sv_type_dict, node_len_dict):

    if a is None:
        return None, None, None

    a_len = 0

    for sv_node in sv_type_dict[sv][0]:

        #Look for aln start node
        if sv_node == aln_nodes[0]:
            a_start = a_len + start - 1
            a_end = a_len + end - 1

        if a == 0:
            #Ref allele : full region length
            a_len += node_len_dict[sv_node]

        else:
            #Alt allele : region length - ref nodes length
            if sv_node not in sv_type_dict[sv][1]:
                a_len += node_len_dict[sv_node]

    return a_len, a_start, a_end

def fill_sv_nodes(ref_path_line):
    pass #implement as function later

def fill_allele_specific_nodes(alt_path_line):
    pass #implement as function later

def check_path_orientation(path, nodes):
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
            return [] ## , ''

    return nodes


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])