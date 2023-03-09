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
REPID = ["A"]
DRAW = [i for i in range(0, 100)][0:10]
DRAW100 = [i for i in range(0, 100)]


rule all:
    input:
        expand(ROOT_dir+"/outputs/simulation/data/{geneID}.phylip", geneID=GENEID),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/pbmpi/{geneID}-{repID}.sh", geneID=GENEID,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/pbmpi/{geneID}-{repID}-cluster.sh", geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/values/{geneID}-{repID}-{draw}.values",geneID=GENEID, repID=REPID, draw=DRAW),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/values/{geneID}-{repID}-{draw}.tre",geneID=GENEID, repID=REPID, draw=DRAW),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.sh", geneID=GENEID, repID=REPID, draw=DRAW),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.fa", geneID=GENEID, repID=REPID, draw=DRAW),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.phylip", geneID=GENEID, repID=REPID, draw=DRAW),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/stats/{geneID}-{repID}-{draw}-OBSERVED.tsv", geneID=GENEID, repID=REPID, draw=DRAW),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/GTR4G/{geneID}-{repID}-{draw}-{repID_}.sh",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/GTR4G/{geneID}-{repID}-{draw}-{repID_}-cluster.sh",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/GTR4G/{geneID}-{repID}-{draw}-{repID_}.sh",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/GTR4G/{geneID}-{repID}-{draw}-{repID_}_ppred.touchdown",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),
        expand(ROOT_dir+"/outputs/simulation/GTR4G/simu/stats/GTR4G/{geneID}-{repID}-{draw}-{repID_}_ppred.tsv",geneID=GENEID, draw=DRAW,repID=REPID,repID_=REPID),

include: ROOT_dir + "/rules/simulation_GTR4G_ppred.smk"
