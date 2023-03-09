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
REPID = ["A"]#["A","B"]
print(GENEID,REPID)

rule all:
    input:
        expand(ROOT_dir+"/outputs/empirical/data/{geneID}.phylip", geneID=GENEID),

        # expand(ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}-cluster.sh", geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.sh", geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.chain", geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}.sh", geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}_ppred.touchdown", geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/ppred_test/stats/{geneID}-GTRG4-{repID}_ppred.tsv",geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/empirical/ppred_test/stats/{geneID}-1_0-OBSERVED.tsv",geneID=GENEID)

include: ROOT_dir + "/rules/empirical_ppred.smk"
