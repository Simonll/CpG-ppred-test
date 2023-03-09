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
GENEID = config["empirical"]["geneID"]
REPID = ["A"]
MODELID = ["GTR4GNT", "F1X4"]

rule all:
    input:
        expand(ROOT_dir+"/outputs/empirical/data/{geneID}.phylip", geneID=GENEID),
        expand(ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}-cluster.sh", geneID=GENEID, modelID=MODELID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}.sh", geneID=GENEID, modelID=MODELID, repID=REPID),
        ## ACTION REQUIRED: to execute pbmpi ##
        expand(ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}.chain", geneID=GENEID, modelID=MODELID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/readpb/{modelID}/{geneID}-{repID}.sh", geneID=GENEID, modelID=MODELID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/readpb/{modelID}/{geneID}-{repID}_ppred.touchdown", geneID=GENEID, modelID=MODELID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/stats/{modelID}/{geneID}-{repID}_ppred.tsv",geneID=GENEID, modelID=MODELID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/stats/{modelID}/{geneID}-1_0-OBSERVED.tsv",geneID=GENEID, modelID=MODELID)

include: ROOT_dir + "/rules/empirical_GTR4G_ppred.smk"
