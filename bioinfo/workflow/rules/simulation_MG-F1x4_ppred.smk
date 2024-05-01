import os
import sys
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)

def genererate_metadata(wildcards):
    geneID = wildcards.geneID
    repID = wildcards.repID
    omega = wildcards.omega
    CpG = wildcards.CpG
    TpA = wildcards.TpA
    tbl = wildcards.tbl
    draw = wildcards.draw
    s:str = '"{'
    s += "'geneID':'"+ geneID+"','repID':'"+repID+"','omega':'"+omega+"','CpGf':'"+CpG+"','TpAf':'"+TpA+"','tbl':'"+tbl+"','draw':'"+draw+"'"+'"}'
    return s

def get_modelID(wildcards)-> str:
    modelID = "MG-F1x4"
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
        docker=ROOT_dir+"/outputs/simulation/MG-F1x4/pbmpi/{geneID}-{repID}.sh",
        cluster=ROOT_dir+"/outputs/simulation/MG-F1x4/pbmpi/{geneID}-{repID}-cluster.sh"
    params:
        local_root=ROOT_dir,
        docker_root="/data",
        cluster_root="/home/sll/CpG-ppred-test/bioinfo/workflow",
        work_dir="/outputs/simulation/MG-F1x4/pbmpi/",
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


rule generate_modified_chain:
    input:
        script=ROOT_dir+"/scripts/generate_MG-F1x4_fixed.py",
        phylip=ROOT_dir+"/outputs/simulation/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/simulation/MG-F1x4/pbmpi/{geneID}-{repID}.chain",
    output:
        ROOT_dir+"/outputs/simulation/MG-F1x4/modified_chain/{geneID}-{omega}-{repID}-{draw}.chain"
    params:
        burnin=100,
        omega="{omega}"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script}\
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --burnin={params.burnin} \
        --omega={params.omega} \
        --output={output}
        """

rule generate_simu_cmd:
    input:
        phylip=ROOT_dir+"/outputs/simulation/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/simulation/MG-F1x4/modified_chain/{geneID}-{omega}-{repID}-{draw}.chain"
    output:
        sh=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}.sh",
        conf=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}.conf"
    params:
        ss="'pwAC pwAG pwAT pwCG pwCT pwGT dinuc31TA dinuc31CA dinuc31TG dinuc31CG nuc3A nuc3C nuc3G nuc3T pwAA'",
        params="'chainID root lambda_CpG lambda_TpA lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT'",
        map="'Nsub Nsynsub dinucCGCA dinucCGTG dinucNSCGCA dinucNSCGTG dinucSCGCA dinucSCGTG gtnrAA gtnrAC gtnrAG gtnrAT gtnrCA gtnrCC gtnrCG gtnrCT gtnrGA gtnrGC gtnrGG gtnrGT gtnrTA gtnrTC gtnrTG gtnrTT'",
        localparams="'-iscodon -code Universal -freeroot Echinops Procavia -rootlength 100 -fixlambdaCpG {CpG} -fixlambdaTpA {TpA} -tofasta'",
        model="M7",
        nsimu="1",
        local=ROOT_dir,
        docker="/data"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 scripts/generate_M0M7GTR_simu_cmd.py \
        --phylip={input.phylip} \
        --mcmc={input.mcmc} \
        --ss={params.ss} \
        --map={params.map} \
        --params={params.params} \
        --localparams={params.localparams} \
        --model={params.model} \
        --nsimu={params.nsimu} \
        --output={output.sh} \
        --local={params.local} \
        --docker={params.docker}
        """

rule run_simu_cmd:
    input:
        ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}}.sh"
    output:
        simu=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}.simu",
        fasta=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.fasta"
    threads: 1
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule fasta2phylip_on_simu:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.fasta"
    output: ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.phylip",
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
        fasta=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.fasta"
    output:
        ROOT_dir+"/outputs/simulation/MG-F1x4/simu/stats/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0-OBSERVED.tsv"
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
        phylip=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.phylip",
    output:
        docker=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}.sh",
        cluster=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}-cluster.sh"
    params:
        local_root=ROOT_dir,
        docker_root="/data",
        cluster_root="/home/sll/CpG-ppred-test/bioinfo/workflow/",
        work_dir="/outputs/simulation/MG-F1x4/simu/pbmpi/",
        phylip="/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        model=get_modelID,
        sampling="'-s -x 10 200'",
        np="4",
        sh_cluster="{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}-cluster.sh",
        sh_docker="{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}.sh",

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
        data=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}.chain"
    output:
        sh=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}.sh",
    params:
        local=ROOT_dir,
        docker="/data",
        sampling="'100 10 200'",
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
        sh=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}.sh",
        phylip=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.phylip",
        param=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}.param",
    output:
        touchdown=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_ppred.touchdown",
        output_2=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_ppred.zip"
    params:
        sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/simulation/MG-F1x4/simu/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-1_0.phylip")+"/'",
        zip_distat_cmd="zip " +  ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_ppred.zip " +  ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_*.ali",
        mv_distat_cmd="mv "+ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_ppred.zip "+ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/",
        rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/simulation/MG-F1x4/simu/pbmpi/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_*.ali",

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
        zip=ROOT_dir+"/outputs/simulation/MG-F1x4/simu/readpb/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_ppred.zip"
    output:
        ROOT_dir+"/outputs/simulation/MG-F1x4/simu/stats/{geneID}-{omega}-{CpG}-{TpA}-{tbl}-{repID}-{draw}-{repID_}_ppred.tsv"
    conda:
        ROOT_dir+"/envs/env.yml"
    params:
        metadata=genererate_metadata
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.zip} --output={output} --metadata={params.metadata}
        """
