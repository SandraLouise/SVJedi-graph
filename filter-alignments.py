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
        which_region = {}
        start_end_region = {}
        rename_node = {}

        with open(gfa_file, "r") as graph_file:
            for line in graph_file:
                if line.startswith("P"):
                    region, nodes, __ = line.split("\t")[1:]
                    chrom = region.split("_")[1]

                    # Add SVs to dict_SVs
                    # if chrom not in dict_SVs.keys():
                    #     dict_SVs[chrom] = []
                    # dict_SVs[chrom].extend(region.split("_")[2:])

                    # Add nodes and region to which_region
                    nodes = extract_nodes(nodes)
                    for node in nodes:
                        which_region[node] = region        

                    # Add first and last node to start_end_region
                    start_end_region[region] = [nodes[0], nodes[-1]]
                
                elif line.startswith("S"):
                    node = line.split("\t")[1]
                    if node[-1] == "a":
                        for n_id in list(which_region.keys()):
                            if n_id.startswith(node[:-1]):
                                which_region[node] = which_region[n_id]
                                break
                            elif n_id.startswith("{}:{}".format(node.split(":")[0], str(int(node.split(":")[1][:-1]) - 1))):
                                new_id = "{}:{}a".format(node.split(":")[0], str(int(node.split(":")[1][:-1]) - 1))
                                which_region[new_id] = which_region[n_id]
                                rename_node[node] = new_id
                                # print(node, n_id, new_id)
                                break
        
        with open(gfa_info_dict, 'w') as file:
            gfa_info = {"which_region" : which_region, "rename_node" : rename_node, "start_end_region" : start_end_region}
            file.write(json.dumps(gfa_info, sort_keys=True, indent=4))
    
    else:
        with open(args.gfainfo[0], "r") as file:
            gfa_info = json.load(file)
            which_region = gfa_info["which_region"]
            start_end_region = gfa_info["start_end_region"]
            rename_node = gfa_info["rename_node"]


    #2. Reading alignments and filling dictionary of informative alignments

    print("Reading gaf file...")

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
            multi_regions = False

            #correct ins nodes id
            for i in range(len(target_nodes)):
                if target_nodes[i] not in which_region.keys():
                    # target_nodes[i] = rename_node[target_nodes[i]]
                    which_region[target_nodes[i]] = which_region[rename_node[target_nodes[i]]]
            
            target_region = which_region[target_nodes[0]] # target_region = "ref"_chrom_ + "_".join([sv_id])
                
            for node in target_nodes[1:]:
                if node not in which_region.keys():
                    # print("Region of {} not found. Path is: {}".format(node, aln["Tid"]))
                    continue
                if which_region[node] != target_region:
                    multi_regions = True
                    
            
            if multi_regions:
                continue
                #TODO: handle multi target-regions (needed for TRANSLOCATIONS)

            else:
                if not is_semiglobal_aln(aln, target_nodes, start_end_region[target_region], d_end):
                    continue
                # print("Passed semiglobal filter:\n", aln)

            # Convert path_based aln to sv-based aln
            target_svs = get_target_svs(target_nodes, target_region)
            # print(target_region)

            for sv_id in target_svs: #sv_id = chrom + "_" + sv_type + "-" + pos + "-" + end
                
                breakpoints = get_breakpoints(sv_id)
                # print(breakpoints)

                #identify supported allele
                allele = get_allele(breakpoints, aln["Tid"], target_nodes)
                # print("Allele: "+str(allele))

                if allele is None:
                    continue

                # BREAKPOINT FILTER
                if not overlap_breakpoints(breakpoints, sv_id, allele, target_nodes, d_over, aln):
                    continue
                
                # print("Passed breakpoint filter")

                # Add aln to dict and save identity
                dict_sv_id = format_dict_sv_id(sv_id)

                if dict_sv_id not in dict_of_informative_aln.keys():
                    dict_of_informative_aln[dict_sv_id] = [[], []]

                dict_of_informative_aln[dict_sv_id][allele].append(line.split("cg:Z:")[0])

                if not identity_saved:
                    list_of_identities.append(aln["Aid"])
                    identity_saved = True
            
    #3. IDENTITY FILTER
    removed_alns = []

    min_id_value = round(stats.quantiles(list_of_identities, n=4)[0], 3)
    # min_id_value = 0.7

    for sv, allele_alns in dict_of_informative_aln.items():
        for allele in [0, 1]:

            curated_alns = []
            for aln in allele_alns[allele]:
                identity = float(aln.rstrip().split("\t")[14].split(":")[-1])

                if identity > min_id_value:
                    curated_alns.append(aln)
                else:
                    removed_alns.append((aln, allele))

            dict_of_informative_aln[sv][allele] = curated_alns
    
    #Stats
    aln_nb = 0
    for sv in dict_of_informative_aln:
        aln_nb += len(dict_of_informative_aln[sv][0]) + len(dict_of_informative_aln[sv][1])

    with open(output_aln_dict, 'w') as file:
        file.write(json.dumps(dict_of_informative_aln, sort_keys=True, indent=4))
    print(str(aln_nb), "informative sv-based aln saved")
    # print(len(removed_alns), " alns removed")



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
        
            for i in range(len(nodes)-1):
           
                # bkpt = [chrom:pos_left, chrom:pos_right]
                chr1, pos1, orient1 = bkpt[0].split(":")
                chr2, pos2, orient2 = bkpt[1].split(":")

                # Check if bkpt link in aln path
                left = all([nodes[i].startswith(chr1), get_node_end(nodes[i]) == pos1, path.split(nodes[i])[0][-1] == orient1])
                right = all([nodes[i+1].startswith(chr2), get_node_start(nodes[i+1]) == pos2, path.split(nodes[i+1])[0][-1] == orient2])

                # Check in bkpt link in aln path (reverse)
                left_rev = all([nodes[i+1].startswith(chr1), get_node_end(nodes[i+1]) == pos1, path.split(nodes[i+1])[0][-1] == invert_orient(orient1)])
                right_rev = all([nodes[i].startswith(chr2), get_node_start(nodes[i]) == pos2, path.split(nodes[i])[0][-1] == invert_orient(orient2)])

                if (left and right) or (left_rev and right_rev):
                    return allele

    
    return None


def format_DEL_INS_id(sv_type, pos, end):
    return '-'.join([sv_type, str(pos), str(end)])
def format_INV_id(pos, end):
    return '-'.join(["INV", str(pos), str(end)])
def format_BND_id(pos, end_chr, end_pos):
    return '-'.join(["BND", str(pos), ":".join([end_chr, end_pos])])

def format_dict_sv_id(graph_sv_id):
    '''
    Graph sv_id:
    - DEL: chrom + "_" + sv_type + "-" + pos + "-" + len
    - INS: chrom + "_" + sv_type + "-" + pos + "-" + len
    - INV: chrom + "_" + sv_type + "-" + pos + "-" + end
    - BND: chrom + "_" + sv_type + "-" + pos + "-" + end_chr + ":" + end_pos

    Dict sv_id:
    - DEL/INS/INV: chrom + "_" + pos + "-" + length
    - BND: ?
    '''
    chrom, info = graph_sv_id.split("_")
    sv_type, pos, end = info.split("-")

    if sv_type in ["DEL", "INS"]:
        return "_".join([chrom, graph_sv_id.split(sv_type + "-")[1]])

    elif sv_type == "INV":
        return "_".join([chrom, "-".join([pos, "0"])])

    # elif type == "BND":
    #     return ""

def get_breakpoints(sv_id):
    '''Returns all breakpoints coordinates (as a list of one duple per bkpt) for a SV (sorted by allele)'''
    #Considering 2 positions possible for left coord (SV start) because of particularities of real SV evaluation set
    
    chrom, info = sv_id.split("_")
    svtype, pos, end = info.split("-")
    bkpts = [[], []] # [[allele 0], [allele 1]]
    forward = ">"
    reverse = "<"

    if svtype == "DEL":
        # end = str(int(pos) + int(end) + 1)
        bkpts[1] = [(":".join([chrom, pos, forward]), ":".join([chrom, end, forward])), 
                    (":".join([chrom, str(int(pos)-1), forward]), ":".join([chrom, str(int(end)-1), forward])), 
                    (":".join([chrom, str(int(pos)-1), forward]), ":".join([chrom, end, forward])),
                    (":".join([chrom, pos, forward]), ":".join([chrom, str(int(end)+1), forward]))] #left coord, right coord
        bkpts[0] = [(":".join([chrom, pos, forward]), ":".join([chrom, str(int(pos)+1), forward])), 
                    (":".join([chrom, str(int(pos)-1), forward]), ":".join([chrom, pos, forward])), 
                    (":".join([chrom, str(int(end)-1), forward]), ":".join([chrom, end, forward])), 
                    (":".join([chrom, str(int(end)-2), forward]), ":".join([chrom, str(int(end)-1), forward]))]
    
    elif svtype == "INS":
        end = str(int(pos) + 1)
        bkpts[0] = [(":".join([chrom, pos, forward]), ":".join([chrom, end, forward])), 
                    (":".join([chrom, str(int(pos)-1), forward]), ":".join([chrom, str(int(end)-1), forward]))] #left coord, right coord
        bkpts[1] = [(":".join([chrom, pos, forward]), ":".join([chrom, str(int(pos)+1)+"a", forward])), 
                    (":".join([chrom, str(int(pos)-1), forward]), ":".join([chrom, str(int(pos)+1)+"a", forward])), 
                    (":".join([chrom, str(int(pos)-1), forward]), ":".join([chrom, pos+"a", forward])), 
                    (":".join([chrom, end+"a", forward]), ":".join([chrom, end, forward])), 
                    (":".join([chrom, end+"a", forward]), ":".join([chrom, str(int(end)-1), forward]))]
    
    elif svtype == "INV":
        bkpts[0] = [(":".join([chrom, pos, forward]), ":".join([chrom, str(int(pos)+1), forward])), 
                    (":".join([chrom, str(int(end)-1), forward]), ":".join([chrom, end, forward]))]
        bkpts[1] = [(":".join([chrom, pos, forward]), ":".join([chrom, str(int(pos)+1), reverse])), 
                    (":".join([chrom, str(int(end)-1), reverse]), ":".join([chrom, end, forward]))]


    elif svtype == "BND":
        pass

    return bkpts

def overlap_breakpoints(bkpts, sv_id, allele, aln_nodes, d_over, aln):

    if any([check_single_breakpoint(sv_id, bkpt, aln_nodes, d_over, aln) for bkpt in bkpts[allele]]):
        return True

    return False

def check_single_breakpoint(sv_id, bkpt, aln_nodes, d_over, aln):
    bkpt_leftCoord = bkpt[0].split(":")[1]
    bkpt_rightCoord = bkpt[1].split(":")[1]
    unaligned_start = aln["Ts"]
    unaligned_end = aln["Tlen"] - aln["Te"] - 1

    #Check presence of breakpoint nodes
    found_leftNode = False
    found_rightNode = False
    for node in aln_nodes:
        if get_node_end(node) == bkpt_leftCoord:
            found_leftNode = True
            bkpt_leftNode = node
        elif get_node_start(node) == bkpt_rightCoord:
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
        # if sum(get_node_len(node) for node in aln_nodes[:aln_nodes.index(bkpt_leftNode)+1]) - unaligned_left >= d_over:
        #     left_overlap = True
        right_overlap = sum(get_node_len(node, sv_id) for node in aln_nodes[aln_nodes.index(bkpt_rightNode):]) - unaligned_end >= d_over
        # if sum(get_node_len(node) for node in aln_nodes[aln_nodes.index(bkpt_rightNode):]) - unaligned_right >= d_over:
        #     right_overlap = True
    
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

def get_target_svs(target_nodes, target_region):
    target_svs = []
    left_node_start = get_node_start(target_nodes[0])
    if left_node_start.endswith("a"):
        left_node_start = int(left_node_start[:-1])
    else:
        left_node_start = int(left_node_start)

    right_node_end = get_node_end(target_nodes[-1])
    if right_node_end.endswith("a"):
        right_node_end = int(right_node_end[:-1])
    else:
        right_node_end = int(right_node_end)

    for sv_id in target_region.split("_")[2:]:
        if not sv_id.startswith("BND"):
            sv_type, sv_start, sv_end = sv_id.split("-")
            # if sv_type in ["DEL", "INS"]:
            #     sv_end = int(sv_start) + int(sv_end)

            if any([min(left_node_start, right_node_end) <= int(sv_start) <= max(left_node_start, right_node_end) , min(left_node_start, right_node_end) <= int(sv_end) <= max(left_node_start, right_node_end)]):
                chrom = target_region.split("_")[1]
                target_svs.append("_".join([chrom, sv_id])) #add chrom to id for later save in aln dict

        else:
            pass
            #TO-DO: handle BNDs

    return target_svs

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
