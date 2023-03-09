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
OMEGA = [0.2,0.8,1]#1
CPG = [1,4,8]
TBL = [1,10]#,2,5,10,
TPA = [1]#4,8
rule all:
    input:
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/data/{geneID}.phylip", geneID=GENEID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.sh", geneID=GENEID,repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr/{geneID}-M0GTR-{repID}-cluster.sh", geneID=GENEID, repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/modified_chain/{geneID}-M0GTR-{omega}-{draw}-{repID}.chain", geneID=GENEID, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh", geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.conf",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.simu",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.fasta",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        expand(ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/stats/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0-OBSERVED.tsv",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-cluster.sh",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.touchdown",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW,repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW, repID=REPID),
        # expand(ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/stats/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.tsv",geneID=GENEID, CpG=CPG, TpA=TPA, tbl=TBL, omega=OMEGA, draw=DRAW, repID=REPID),

include: ROOT_dir + "/rules/simulation_m0gtr_ppred.smk"
