import os
import sys
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)


def genererate_metadata(wildcards):
    geneID = wildcards.geneID
    repID = wildcards.repID
    draw = wildcards.draw
    s:str = '"{'
    s += "'geneID':'"+ geneID+"','repID':'"+repID+"','draw':'"+draw+"'"+'"}'
    return s

def get_modelID(wildcards)-> str:
    modelID = "GTR4GNT"
    print(config["modelID"][modelID])
    return "'"+config["modelID"][modelID]+"'"

def get_Nsite(wildcards)-> int:
    geneID = wildcards.geneID
    return int(config["Nsites"][geneID])

rule fasta2phylip:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/data/{geneID}.fasta"
    output:
        ROOT_dir+"/outputs/simulation/data/{geneID}.phylip"
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
        phylip=ROOT_dir+"/outputs/simulation/data/{geneID}.phylip",
    output:
        docker=ROOT_dir+"/outputs/simulation/GTR4G/pbmpi/{geneID}-{repID}.sh",
        cluster=ROOT_dir+"/outputs/simulation/GTR4G/pbmpi/{geneID}-{repID}-cluster.sh"
    params:
        local_root=ROOT_dir,
        docker_root="/data",
        cluster_root="/home/sll/CpG-ppred-test/bioinfo/workflow",
        work_dir="/outputs/simulation/GTR4G/pbmpi/",
        phylip="/outputs/simulation/data/{geneID}.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
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
        --sh_docker={params.sh_docker}
        """

rule generate_GTR4G_values:
    input:
        script=ROOT_dir+"/scripts/generate_GTR4G_values.py",
        mcmc=ROOT_dir+"/outputs/simulation/GTR4G/pbmpi/{geneID}-{repID}.chain",
    output:
        values=ROOT_dir+"/outputs/simulation/GTR4G/values/{geneID}-{repID}-{draw}.values",
        tree=ROOT_dir+"/outputs/simulation/GTR4G/values/{geneID}-{repID}-{draw}.tre"
    params:
        burnin=100,
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script}\
        --mcmc={input.mcmc} \
        --burnin={params.burnin} \
        --output_values={output.values} \
        --output_tree={output.tree}
        """

rule generate_simu_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_GTR4G_simu_cmd.py",
        values=ROOT_dir+"/outputs/simulation/GTR4G/values/{geneID}-{repID}-{draw}.values",
        tree=ROOT_dir+"/outputs/simulation/GTR4G/values/{geneID}-{repID}-{draw}.tre"
    output:
        sh=ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.sh"
    params:
        simu=ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}",
        length=get_Nsite,
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script}\
        --simu={params.simu} \
        --values={input.values} \
        --tree={input.tree} \
        --length={params.length} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_simu_cmd:
    input:
        ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.sh"
    output:
        simu=ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.fa",
    threads: 1
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule fasta2phylip_simu:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.fa"
    output:
        ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.phylip"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.data} \
        --phylip={output}
        """


rule compute_CpGfreq_realdata:
    input:
        script=ROOT_dir+"/scripts/compute_CpGfreq.py",
        fasta=ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.fa"
    output:
        ROOT_dir+"/outputs/simulation/GTR4G/simu/stats/{geneID}-{repID}-{draw}-OBSERVED.tsv"
    conda:
        ROOT_dir+"/envs/env.yml"

    params:
        metadata=genererate_metadata
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.fasta} --output={output} --metadata={params.metadata}
        """

rule generate_pb_mpi_on_simu_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.phylip",
    output:
        docker=ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}.sh",
        cluster=ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}-cluster.sh"
    params:
        local_root=ROOT_dir,
        docker_root="/data",
        cluster_root="/home/sll/CpG-ppred-test/bioinfo/workflow/",
        work_dir="/outputs/simulation/GTR4G/simu/pbmpi/",
        phylip="/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        model=get_modelID,
        sampling="'-s -x 10 200'",
        np="4",
        sh_cluster="{geneID}-{repID}-{draw}-cluster.sh",
        sh_docker="{geneID}-{repID}-{draw}.sh",

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
        --sh_docker={params.sh_docker}
        """

rule generate_readpb_mpi_ppred_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_readpb_mpi_ppred_cmd.py",
        data=ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}.chain"
    output:
        sh=ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/{geneID}-{repID}-{draw}-{repID_}.sh",
    params:
        local=ROOT_dir,
        docker="/data",
        sampling="'100 10 200'",
        model="GTRG4",
        image="'ubuntu20.04/pbmpi_mapstats:latest'",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --input={input.data} \
        --output={output.sh} \
        --model={params.model} \
        --image={params.image} \
        --local={params.local} \
        --docker={params.docker} \
        --sampling={params.sampling} \
        """

rule run_readpb_mpi_ppred_cmd:
    input:
        sh=ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/{geneID}-{repID}-{draw}-{repID_}.sh",
        phylip=ROOT_dir+"/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.phylip",
        param=ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}-{repID_}.param",
    output:
        touchdown=ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/{geneID}-{repID}-{draw}-{repID_}_ppred.touchdown",
        output_2=ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/{geneID}-{repID}-{draw}-{repID_}_ppred.zip"
    params:
        sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/simulation/GTR4G/simu/{geneID}-{repID}-{draw}.phylip")+"/'",
        zip_distat_cmd="zip " +  ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}-{repID_}_ppred.zip " +  ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}-{repID_}_*.ali",
        mv_distat_cmd="mv "+ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}-{repID_}_ppred.zip "+ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/",
        rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/simulation/GTR4G/simu/pbmpi/{geneID}-{repID}-{draw}-{repID_}_*.ali",

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
        zip=ROOT_dir+"/outputs/simulation/GTR4G/simu/readpb/{geneID}-{repID}-{draw}-{repID_}_ppred.zip"
    output:
        ROOT_dir+"/outputs/simulation/GTR4G/simu/stats/{geneID}-{repID}-{draw}-{repID_}_ppred.tsv"
    conda:
        ROOT_dir+"/envs/env.yml"
    params:
        metadata=genererate_metadata
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.zip} --output={output} --metadata={params.metadata}
        """
