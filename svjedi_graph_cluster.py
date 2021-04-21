#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

def get_script_path():
    '''source: https://stackoverflow.com/a/4943474'''
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def main(args):

    parser = argparse.ArgumentParser(
        description="---"
    )    
    parser.add_argument(
        "-v", 
        "--vcf", 
        metavar="<inputVCF_absolutePath>", 
        type=str,
        nargs=1, 
        required=True)

    parser.add_argument(
        "-r", 
        "--ref", 
        metavar="<inputReferenceGenome_absolutePath>", 
        type=str,
        nargs=1, 
        required=True)
    
    parser.add_argument(
        "-f", 
        "--readfile", 
        metavar="<inputReadFile_absolutePath>", 
        type=str,
        nargs=1, 
        required=True)
    
    parser.add_argument(
        "-V", 
        "--version", 
        metavar="<scriptVersionToUse", 
        type=str,
        required=False)
    
    parser.add_argument(
        "-o", 
        "--output", 
        metavar="<output>", 
        type=str,
        required=False)

    args = parser.parse_args()

    inputVCF = args.vcf[0]
    inputREF = args.ref[0]
    inputREAD = args.readfile[0]
    
    if args.version and args.version[0] in ["3", "3.5"]:
        data_prep_script = "{}/extract_variant_data_v3.py".format(get_script_path())
        aln_filter_script = "{}/filter_aln_v{}.py".format(get_script_path(), args.version[0])
        version = "v" + args.version[0]
    else:
        data_prep_script = "{}/extract_variant_data_v3.py".format(get_script_path())
        aln_filter_script = "{}/filter_aln_v3.py".format(get_script_path())
        version = "v3"
    geno_pred_script = "{}/svjedi_genotype_prediction.py".format(get_script_path())
    
    if args.output:
        output_label = args.output[0]
    else:
        #unique label to avoid unintentional overwriting of previous results
        import time
        tpt = {}
        tpt["yr"], tpt["mth"], tpt["dy"], tpt["hr"], tpt["mn"], tpt["sc"], *__ = list(map(str, time.localtime()))
        for key, val in tpt.items():
            if len(val) == 1:
                tpt[key] = "0" + val
        output_label = "-".join(["".join([tpt["yr"], tpt["mth"], tpt["dy"]]), "".join([tpt["hr"], tpt["mn"], tpt["sc"]])])
    

    ## SVJedi-Graph steps
    cmds = []

    #0. Create dir for intermediate files and final results
    dir_name = "svj-g-{}_{}".format(version, output_label)
    cmds.append("mkdir {}".format(dir_name))
    cmds.append("cd {}".format(dir_name))

    #1. Extract variant data
    #load python env
    cmds.append(". /softs/local/env/envpython-3.8.5.sh")
    #
    cmds.append("python3 {} -v {} -r {} -o {}".format(data_prep_script, inputVCF, inputREF, dir_name))
    outputVCF = "{}/converted_variants.vcf".format(dir_name)
    outputREF = "{}/reference_subsequences.fa".format(dir_name)

    #2. Create graph
    #load conda env
    cmds.append(". /local/env/envconda.sh")
    #
    outputGRAPH = "{}/variant_graph".format(dir_name)
    cmds.append("conda activate env_vg")
    cmds.append("~/.conda/envs/env_vg/bin/vg construct -r {} -v {} -m 5000 -S -a -f > {}.vg".format(outputREF, outputVCF, outputGRAPH))
    cmds.append("~/.conda/envs/env_vg/bin/vg view --gfa {}.vg > {}.gfa".format(outputGRAPH, outputGRAPH))
    cmds.append("conda deactivate")

    #3. Map reads
    #graphaligner version
    graphaligner_env = "~/.conda/envs/graphaligner_1.0.13"
    #
    outputGAF = "{}/raw_aln.gaf".format(dir_name)
    cmds.append("conda activate {}".format(graphaligner_env))
    cmds.append("{}/bin/GraphAligner -g {}.gfa -f {} -a {} -x vg > {}/mapping_log.txt".format(graphaligner_env, outputGRAPH, inputREAD, outputGAF, dir_name))
    cmds.append("conda deactivate")

    #4. Filter alns
    #load python env
    cmds.append(". /softs/local/env/envpython-3.8.5.sh")
    #
    cmds.append("python3 {} -a {} -g {}.gfa".format(aln_filter_script, outputGAF, outputGRAPH))
    outputJSON = "{}/informative_aln_{}.json".format(dir_name, version)

    #5. Predict genotype
    cmds.append("python3 {} -d {} -v {}".format(geno_pred_script, outputJSON, inputVCF))

    ## Running SVJedi-Graph
    for cmd in cmds:
        subprocess.run(cmd, shell=True)
    

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])