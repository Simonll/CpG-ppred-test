import os
import sys
import re
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)


def get_modelID(wildcards)-> str:
    modelID = wildcards.modelID
    print(config["modelID"][modelID])
    return "'"+config["modelID"][modelID]+"'"

def genererate_metadata(wildcards):
    s:str = '"{'
    geneID = wildcards.geneID
    modelID = wildcards.modelID
    s += "'geneID':'"+ geneID+"'"
    s += ",'modelID':'"+modelID+"'"
    if "repID" in wildcards:
        repID = wildcards.repID
        s+= ",'repID':'"+repID+"'"
    s+='"}'
    return s

def get_tree(wildcards):
    if wildcards.geneID == "concat":
        return ROOT_dir+"/data/"+"Mammalia-39sp-CAT_unrooted_withoutsupport.tre"

def get_tree_docker(wildcards):
    if wildcards.geneID == "concat":
        return "/data/"+"Mammalia-39sp-CAT_unrooted_withoutsupport.tre"


rule fasta2phylip:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/data/{geneID}.fasta"
    output: ROOT_dir+"/outputs/empirical/data/{geneID}.phylip"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.data} \
        --phylip={output}
        """

rule generate_pb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/empirical/data/{geneID}.phylip",
    output:
        docker=ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}.sh",
        cluster=ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}-cluster.sh"
    params:
        local_root=ROOT_dir,
        docker_root="/data",
        cluster_root="/home/sll/CpG-ppred-test/bioinformatics/workflow/",
        work_dir="/outputs/empirical/pbmpi/{modelID}/",
        phylip="/outputs/empirical/data/{geneID}.phylip",
        tree=get_tree,
        model=get_modelID,
        sampling="'-s -x 10 200'",
        np="4",
        sh_cluster="{geneID}-{repID}-cluster.sh",
        sh_docker="{geneID}-{repID}.sh",

    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={params.phylip} \
        --tree={params.tree} \
        --model={params.model} \
        --sampling={params.sampling} \
        --np={params.np} \
        --local_root={params.local_root} \
        --docker_root={params.docker_root} \
        --cluster_root={params.cluster_root} \
        --work_dir={params.work_dir} \
        --sh_cluster={params.sh_cluster} \
        --sh_docker={params.sh_docker} \
        """

rule generate_readpb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_readpb_mpi_ppred_cmd.py",
        data=ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/empirical/readpb/{modelID}/{geneID}-{repID}.sh",
    params:
        local=ROOT_dir,
        docker="/data",
        sampling="'100 1 200'",
        image="'ubuntu20.04/pbmpi:latest'",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input={input.data} \
        --output={output.sh} \
        --image={params.image} \
        --local={params.local} \
        --docker={params.docker} \
        --sampling={params.sampling} \
        """

rule run_readpb_mpi_cmd:
    input:
        sh=ROOT_dir+"/outputs/empirical/readpb/{modelID}/{geneID}-{repID}.sh",
        phylip=ROOT_dir+"/outputs/empirical/data/{geneID}.phylip",
        param=ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}.param",
    output:
        touchdown=ROOT_dir+"/outputs/empirical/readpb/{modelID}/{geneID}-{repID}_ppred.touchdown",
        output_2=ROOT_dir+"/outputs/empirical/readpb/{modelID}/{geneID}-{repID}_ppred.zip"
    params:
        sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/empirical/data/{geneID}.phylip")+"/'",
        zip_distat_cmd="zip " +  ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}_ppred.zip " +  ROOT_dir+"/outputs/empirical/pbmpi//{modelID}/{geneID}-{repID}_ppred*.ali",
        mv_distat_cmd="mv "+ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}_ppred.zip "+ROOT_dir+"/outputs/empirical/readpb/{modelID}/",
        rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/empirical/pbmpi/{modelID}/{geneID}-{repID}_*.ali",
    conda:
        ROOT_dir+"/envs/env.yml"
    threads: 2
    shell:
        """
        {params.sed_cmd} {input.param} && bash {input.sh} && \
        {params.zip_distat_cmd} && {params.rm_distat_s_cmd} && {params.mv_distat_cmd} && \
        touch {output.touchdown}
        """

rule compute_CpGfreq:
    input:
        script=ROOT_dir+"/scripts/compute_CpGfreq.py",
        zip=ROOT_dir+"/outputs/empirical/readpb/{modelID}/{geneID}-{repID}_ppred.zip"
    output:
        ROOT_dir+"/outputs/empirical/stats/{modelID}/{geneID}-{repID}_ppred.tsv"
    params:
        metadata=genererate_metadata,
    conda:
        ROOT_dir+"/envs/env.yml"
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.zip} --output={output} --metadata={params.metadata}
        """

rule compute_CpGfreq_realdata:
    input:
        script=ROOT_dir+"/scripts/compute_CpGfreq.py",
        fasta=ROOT_dir+"/data/{geneID}.fasta"
    output:
        ROOT_dir+"/outputs/empirical/stats/{geneID}-1_0-OBSERVED.tsv"
    params:
        metadata=genererate_metadata,
    conda:
        ROOT_dir+"/envs/env.yml"
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.fasta} --output={output} --metadata={params.metadata}
        """
