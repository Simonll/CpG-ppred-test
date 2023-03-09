import os
import sys
import re
from pathlib import Path
import json

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).__str__()
workdir: ROOT_dir
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

configfile: ROOT_dir + "/configs/configs.yaml"
GENEID = config["simulation"]["geneID"]
REPID = ["A"]#["A","B"]
DRAW = [i for i in range(0, 100)][0:10]
DRAW100 = [i for i in range(0, 100)]

# ruleorder: fasta2phylip_simu > generate_GTRG4_pb_mpi_cmd > serialize_mapping_pred > generate_two_sites_mappings_pred > generate_two_sites_mapping_stats_pred > aggregate_two_sites_results_pred

rule all:
    input:
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/data/{geneID}.phylip", geneID=GENEID),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.sh", geneID=GENEID,repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4/{geneID}-GTRG4-{repID}-cluster.sh", geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/values/{geneID}-GTRG4-{repID}-{draw}.values",geneID=GENEID, repID=REPID, draw=DRAW),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/values/{geneID}-GTRG4-{repID}-{draw}.tre",geneID=GENEID, repID=REPID, draw=DRAW),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred/values/{geneID}-GTRG4-{repID}-{draw}.values",geneID=GENEID, repID=REPID, draw=DRAW100),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred/values/{geneID}-GTRG4-{repID}-{draw}.tre",geneID=GENEID, repID=REPID, draw=DRAW100),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.sh", geneID=GENEID, repID=REPID, draw=DRAW),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.fa", geneID=GENEID, repID=REPID, draw=DRAW),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred/simu/{geneID}-GTRG4-{repID}-{draw}.sh", geneID=GENEID, repID=REPID, draw=DRAW100),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred/simu/{geneID}-GTRG4-{repID}-{draw}.fa", geneID=GENEID, repID=REPID, draw=DRAW100),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip", geneID=GENEID, repID=REPID, draw=DRAW),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.sh",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}-cluster.sh",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.sh",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_map.touchdown",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        # expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred/stats/{geneID}-GTRG4-{repID}-{draw}.tsv",geneID=GENEID, repID=REPID, draw=DRAW100),
        expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/stats/{geneID}-GTRG4-{repID}-{draw}-OBSERVED.tsv", geneID=GENEID, repID=REPID, draw=DRAW),
        expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.sh",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.touchdown",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        expand(ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/stats/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.tsv",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        # ROOT_di)r+"/outputs/simulation/gtrg4/stats_gtrg4/aggregated_stats_DINT.tsv"
include: ROOT_dir + "/rules/simulation_gtrg4_ppred.smk"
