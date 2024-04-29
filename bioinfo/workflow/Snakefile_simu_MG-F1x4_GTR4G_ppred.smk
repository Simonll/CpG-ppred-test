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
OMEGA = [0.2,1]
CPG = [1,4,8]
TBL = [1,10]
TPA = [1]
rule all:
    input:
        expand(ROOT_dir+"/outputs/simulation/data/{geneID}.phylip", geneID=GENEID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/pbmpi/{geneID}-{repID}.sh", geneID=GENEID,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/pbmpi/{geneID}-{repID}-cluster.sh", geneID=GENEID, repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/modified_chain/{geneID}-{repID}-{omega}-{draw}.chain", geneID=GENEID, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}.sh", geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}.conf",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}.simu",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-1_0.fasta",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-1_0.phylip",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/stats/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-1_0-OBSERVED.tsv",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/GTR4G/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID_}.sh",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/GTR4G/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID_}-cluster.sh",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/GTR4G/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID_}.sh",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/GTR4G/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID_}_ppred.touchdown",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/GTR4G/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID_}_ppred.zip",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW, repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/MG-F1x4/simu/stats/GTR4G/{geneID}-{repID}-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID_}_ppred.tsv",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW, repID=REPID),

include: ROOT_dir + "/rules/simulation_MG-F1x4_ppred.smk"
