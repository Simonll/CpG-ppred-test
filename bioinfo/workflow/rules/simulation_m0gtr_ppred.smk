import os
import sys
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)



def genererate_metadata(wildcards):
    s:str = '"{'
    geneID = wildcards.geneID
    repID = wildcards.repID
    omega = wildcards.omega
    CpG = wildcards.CpG
    TpA = wildcards.TpA
    tbl = wildcards.tbl
    draw = wildcards.draw
    s += "'geneID':'"+ geneID+"','repID':'"+repID+"','omega':'"+omega+"','CpGf':'"+CpG+"','TpAf':'"+TpA+"','tbl':'"+tbl+"','draw':'"+draw+"'"+'"}'
    return s


rule fasta2phylip:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/data/{geneID}.fasta"
    output:
        ROOT_dir+"/outputs/simulation/m0gtr/data/{geneID}.phylip"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.data} \
        --phylip={output}
        """

rule generate_M0GTR_pb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_M0GTR_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/simulation/m0gtr/data/{geneID}.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:
        docker=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.sh",
        cluster=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr/{geneID}-M0GTR-{repID}-cluster.sh"
    params:
        local=ROOT_dir,
        docker="/data",
        phylip="/outputs/simulation/m0gtr/data/{geneID}.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        cluster="/home/sll/mappings_brief",
        chainname="{geneID}-M0GTR-{repID}",
        np="4",
        output="{geneID}-M0GTR-{repID}-cluster.sh"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --tree={input.tree} \
        --output={output.docker} \
        --local={params.local} \
        --docker={params.docker}
        echo #####
        echo "#qsub -q all.q@compute-0-8.local -q all.q@compute-0-9.local -q all.q@compute-0-10.local -q all.q@compute-0-11.local -q all.q@compute-0-12.local -q all.q@compute-0-13.local -q all.q@compute-0-14.local -q all.q@compute-0-15.local -cwd -pe mpi {params.np} {params.output}" >> {output.cluster}
        echo "LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:\$LD_LIBRARY_PATH" >> {output.cluster}
        echo "CXX=\`which g++\`" >> {output.cluster}
        echo "CC=\`which gcc\`" >> {output.cluster}
        echo "export LD_LIBRARY_PATH" >> {output.cluster}
        echo "export CC" >> {output.cluster}
        echo "export CXX" >> {output.cluster}
        echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} -mutsel -catfix uniform -s -x 1 2000 {params.chainname} >> {output.cluster}
        """

rule generate_M0GTR_pb_mpi_cmd_:
    input:
        script=ROOT_dir+"/scripts/generate_M0GTRW_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}-1_0.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:
        docker=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}.sh",
        cluster=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}-cluster.sh"
    params:
        local=ROOT_dir,
        docker="/data",
        phylip="/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}-1_0.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        cluster="/home/sll/mappings_brief",
        chainname="{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}",
        np="4",
        output="{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}-cluster.sh"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --tree={input.tree} \
        --output={output.docker} \
        --local={params.local} \
        --docker={params.docker}
        echo #####
        echo "#qsub -q all.q@compute-0-8.local -q all.q@compute-0-9.local -q all.q@compute-0-10.local -q all.q@compute-0-11.local -q all.q@compute-0-12.local -q all.q@compute-0-13.local -q all.q@compute-0-14.local -q all.q@compute-0-15.local -cwd -pe mpi {params.np} {params.output}" >> {output.cluster}
        echo "LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:\$LD_LIBRARY_PATH" >> {output.cluster}
        echo "CXX=\`which g++\`" >> {output.cluster}
        echo "CC=\`which gcc\`" >> {output.cluster}
        echo "export LD_LIBRARY_PATH" >> {output.cluster}
        echo "export CC" >> {output.cluster}
        echo "export CXX" >> {output.cluster}
        echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} -mutsel -catfix uniform -freeomega -s -x 10 200 {params.chainname} >> {output.cluster}
        """

rule generate_modified_chain:
    input:
        script=ROOT_dir+"/scripts/generate_M0GTR_fixed.py",
        phylip=ROOT_dir+"/outputs/simulation/m0gtr/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.chain",
    output:
        ROOT_dir+"/outputs/simulation/m0gtr/modified_chain/{geneID}-M0GTR-{omega}-{draw}-{repID}.chain"
    params:
        burnin=1000,
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



# rule generate_simu_cmd:
#     input:
#         phylip=ROOT_dir+"/outputs/simulation/m0gtr/data/{geneID}.phylip",
#         mcmc=ROOT_dir+"/outputs/simulation/m0gtr/modified_chain/{geneID}-M0GTR-{omega}-{draw}-{repID}.chain"
#     output:
#         sh=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{draw}-{repID}.sh",
#         conf=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{draw}-{repID}.conf"
#     params:
#         ss="'pwAC pwAG pwAT pwCG pwCT pwGT dinuc31CG dinuc31TG dinuc31CA nuc3A nuc3C nuc3G nuc3T pwAA'",
#         params="'chainID root lambda_CpG lambda_TBL lambda_omega nucsA nucsC nucsG nucsT nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT'",
#         map="'Nsub Nsynsub dinucCGCA dinucCGTG dinucNSCGCA dinucNSCGTG dinucSCGCA dinucSCGTG gtnrAA gtnrAC gtnrAG gtnrAT gtnrCA gtnrCC gtnrCG gtnrCT gtnrGA gtnrGC gtnrGG gtnrGT gtnrTA gtnrTC gtnrTG gtnrTT'",
#         localparams="'-iscodon -code Universal -freeroot Echinops Procavia -rootlength 100 -fixlambdaCpG {CpG} -tofasta'",
#         model="M7",
#         nsimu="1",
#         local=ROOT_dir,
#         docker="/data"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 scripts/generate_M0M7GTR_simu_cmd.py \
#         --phylip={input.phylip} \
#         --mcmc={input.mcmc} \
#         --ss={params.ss} \
#         --map={params.map} \
#         --params={params.params} \
#         --localparams={params.localparams} \
#         --model={params.model} \
#         --nsimu={params.nsimu} \
#         --output={output.sh} \
#         --local={params.local} \
#         --docker={params.docker}
#         """

rule generate_simu_cmd_:
    input:
        phylip=ROOT_dir+"/outputs/simulation/m0gtr/data/{geneID}.phylip",
        mcmc=ROOT_dir+"/outputs/simulation/m0gtr/modified_chain/{geneID}-M0GTR-{omega}-{draw}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}.sh",
        conf=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{draw}-{tbl}-{repID}.conf"
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
        ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh"
    output:
        simu=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.simu",
        fasta=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.fasta"
    threads: 1
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule fasta2phylip_on_simu:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.fasta"
    output: ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip",
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
        fasta=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.fasta"
    output:
        ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/stats/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0-OBSERVED.tsv"
    params:
        metadata=genererate_metadata,
    threads: 1
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} --input={input.fasta} --output={output} --metadata={params.metadata}
        """



rule generate_GTRG4_pb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_GTRG4_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:
        docker=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh",
        cluster=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-cluster.sh"
    params:
        local=ROOT_dir,
        docker="/data",
        phylip="/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        cluster="/home/sll/mappings_brief",
        chainname="{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}",
        np="4",
        output="{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-cluster.sh"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --phylip={input.phylip} \
        --tree={input.tree} \
        --output={output.docker} \
        --local={params.local} \
        --docker={params.docker}
        echo #####
        echo "#qsub -q all.q@compute-0-8.local -q all.q@compute-0-9.local -q all.q@compute-0-10.local -q all.q@compute-0-11.local -q all.q@compute-0-12.local -q all.q@compute-0-13.local -q all.q@compute-0-14.local -q all.q@compute-0-15.local -cwd -pe mpi {params.np} {params.output}" >> {output.cluster}
        echo "LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:\$LD_LIBRARY_PATH" >> {output.cluster}
        echo "CXX=\`which g++\`" >> {output.cluster}
        echo "CC=\`which gcc\`" >> {output.cluster}
        echo "export LD_LIBRARY_PATH" >> {output.cluster}
        echo "export CC" >> {output.cluster}
        echo "export CXX" >> {output.cluster}
        echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} -gtr -ncat 1 -dgam 4 -s -x 10 200 {params.chainname} >> {output.cluster}
        """


rule generate_readpb_mpi_ppred_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_readpb_mpi_ppred_cmd.py",
        data=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh",
    params:
        local=ROOT_dir,
        docker="/data",
        sampling="'100 2 200'",
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

rule generate_readpb_mpi_ppred_cmd_:
    input:
        script=ROOT_dir+"/scripts/generate_readpb_mpi_ppred_cmd.py",
        data=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh",
    params:
        local=ROOT_dir,
        docker="/data",
        sampling="'100 2 200'",
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


rule run_NT_readpb_mpi_ppred_cmd:
    input:
        sh=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh",
        phylip=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip",
        param=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.param",
    output:
        touchdown=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.touchdown",
        output_2=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip"
    params:
        sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip")+"/'",
        zip_distat_cmd="zip " +  ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip " +  ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_*.ali",
        mv_distat_cmd="mv "+ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip "+ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_gtrg4/",
        rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_*.ali",

    threads: 2
    shell:
        """
        {params.sed_cmd} {input.param} && bash {input.sh} && \
        {params.zip_distat_cmd} && {params.rm_distat_s_cmd} && {params.mv_distat_cmd} && \
        touch {output.touchdown}
        """

rule run_NT_readpb_mpi_ppred_cmd_:
    input:
        sh=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.sh",
        phylip=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip",
        param=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}.param",
    output:
        touchdown=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.touchdown",
        output_2=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip"
    params:
        sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}-1_0.phylip")+"/'",
        zip_distat_cmd="zip " +  ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip " +  ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_*.ali",
        mv_distat_cmd="mv "+ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip "+ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/",
        rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_*.ali",

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
        zip=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip"
    output:
        ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/stats/{geneID}-GTRG4-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.tsv"
    params:
        metadata=genererate_metadata,
    threads: 1
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} --input={input.zip} --output={output} --metadata={params.metadata}
        """

rule compute_CpGfreq_:
    input:
        script=ROOT_dir+"/scripts/compute_CpGfreq.py",
        zip=ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/readpb_m0gtr_/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.zip"
    output:
        ROOT_dir+"/outputs/simulation/m0gtr/ppred_test/stats/{geneID}-M0GTR-{omega}-{CpG}-{TpA}-{tbl}-{draw}-{repID}_ppred.tsv"
    params:
        metadata=genererate_metadata,
    threads: 1
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} --input={input.zip} --output={output} --metadata={params.metadata}
        """



# rule compute_stats_DINT:
#     input:
#         script=ROOT_dir+"/scripts/compute_CpG_stat_.py",
#         zipf=ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_suffdistatmap.zip"
#     output:
#         ROOT_dir+"/outputs/simulation/m0gtr/stats_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_suffdistatmap_DINT.tsv",
#     params:
#         metadata=genererate_metadata
#     threads: 1
#     shell:
#         """
#         python3 {input.script} --input={input.zipf} --output={output} --metadata={params.metadata} && touch {output}
#         """

# rule aggregated_stats_DINT:
#     input:
#         script=ROOT_dir+"/scripts/aggregate_CpG_stats.py",
#     output:
#          ROOT_dir+"/outputs/simulation/m0gtr/stats_gtrg4/aggregated_stats_DINT.tsv"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/stats_gtrg4/",
#         pattern="*_suffdistatmap_DINT.tsv",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} --input_dir={params.input_dir} --pattern={params.pattern} --output={output}
#         """

# rule run_NT_readpb_mpi_cmd:
#     input:
#         sh=ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}.sh",
#         phylip=ROOT_dir+"/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{draw}-{repID}-1_0.phylip",
#         param=ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}.param",
#     output:
#         ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_map.touchdown"
#     params:
#         sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/simulation/m0gtr/simu/{geneID}-M0GTR-{omega}-{CpG}-{draw}-{repID}-1_0.phylip")+"/'",
#         rm_cmd_nodestates=get_nodestates_to_delete,
#         mv_cmd="mv "+ROOT_dir+"/outputs/simulation/m0gtr/pbmpi_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_*.map "+ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/",
#     threads: 2
#     shell:
#         """
#         {params.sed_cmd} {input.param} && bash {input.sh}  && touch {output}
#         """
#         #&& {params.rm_cmd_nodestates} && {params.mv_cmd}

# rule serialize_mapping:
#     input:
#         script=ROOT_dir+"/scripts/serialize_mapping_.py",
#     output:
#         ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_*.map'",
#         output_dir=ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/",
#         metadata=get_params,
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --output_dir={params.output_dir} \
#         --pattern={params.pattern} \
#         --metadata={params.metadata} && touch {output}
#         """
#         #&& rm $(echo {params.input_dir}{params.pattern})

# rule generate_two_sites_mappings:
#     input:
#         script=ROOT_dir+"/scripts/generate_two_sites_mapping_.py",
#         prevtouchdown=ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites_/",
#         list_of_paired_sites=get_adjsites
#     output:
#         #temp(expand(ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites_/{{geneID}}-GTRG4-{{omega}}-{{CpG}}-{{draw}}-{{repID}}_{adj_sites}.pickle", adj_sites=ADJ_SITES))
#         ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites_/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites.touchdown"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --list_of_paired_sites={params.list_of_paired_sites} \
#         --output_dir={params.output_dir}  &&  touch {output}
#         """

# rule generate_two_sites_mapping_stats:
#     input:
#         script=ROOT_dir+"/scripts/generate_two_sites_mapping_stats_.py",
#         prevtouchdown=ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites_/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites_/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites/",
#     output:
#         #temp(expand(ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites/{{geneID}}-GTRG4-{{omega}}-{{CpG}}-{{draw}}-{{repID}}_{adj_sites}-{type}.pickle", adj_sites=ADJ_SITES, type=["SUB","TIME"])),
#         sub=ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites-SUB.touchdown",
#         time=ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites-TIME.touchdown",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --output_dir={params.output_dir}  && touch {output.sub} && touch {output.time}
#         """
#         #&& rm $(echo {params.input_dir}{params.pattern})

# rule aggregate_two_sites_results:
#     input:
#         script=ROOT_dir+"/scripts/aggregate_results_over_sites.py",
#         prevtouchdown=expand(ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites/{{geneID}}-GTRG4-{{omega}}-{{CpG}}-{{draw}}-{{repID}}_two_sites-{type}.touchdown", type=["SUB","TIME"])
#     output:
#         sub=ROOT_dir+"/outputs/simulation/m0gtr/stats/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}-two_sites-SUB.pickle",
#         time=ROOT_dir+"/outputs/simulation/m0gtr/stats/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}-two_sites-TIME.pickle",
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/step_two_sites/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}*'"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --output_sub={output.sub} \
#         --output_time={output.time}
#         """
#         #&& rm $(echo {params.input_dir}{params.pattern})

# rule serialize_mapping_pred:
#     input:
#         script=ROOT_dir+"/scripts/serialize_mapping_.py",
#     output:
#         ROOT_dir+"/outputs/simulation/m0gtr/pred/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_*.map'",
#         output_dir=ROOT_dir+"/outputs/simulation/m0gtr/pred/readpb_gtrg4/",
#         metadata=get_params,
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --output_dir={params.output_dir} \
#         --pattern={params.pattern} \
#         --metadata={params.metadata} && touch {output}
#         """
#         #&& rm $(echo {params.input_dir}{params.pattern})
# rule generate_two_sites_mappings_pred:
#     input:
#         script=ROOT_dir+"/scripts/generate_two_sites_mapping_.py",
#         prevtouchdown=ROOT_dir+"/outputs/simulation/m0gtr/pred/readpb_gtrg4/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/pred/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites_/",
#         list_of_paired_sites=get_adjsites
#     output:
#         ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites_/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites.touchdown"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --list_of_paired_sites={params.list_of_paired_sites} \
#         --output_dir={params.output_dir}  &&  touch {output}
#         """

# rule generate_two_sites_mapping_stats_pred:
#     input:
#         script=ROOT_dir+"/scripts/generate_two_sites_mapping_stats_.py",
#         prevtouchdown=ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites_/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites_/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites/",
#     output:
#         sub=ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites-SUB.touchdown",
#         time=ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}_two_sites-TIME.touchdown",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --output_dir={params.output_dir}  && touch {output.sub} && touch {output.time}
#         """
#         #&& rm $(echo {params.input_dir}{params.pattern})


# rule aggregate_two_sites_results_pred:
#     input:
#         script=ROOT_dir+"/scripts/aggregate_results_over_sites.py",
#         prevtouchdown=expand(ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites/{{geneID}}-GTRG4-{{omega}}-{{CpG}}-{{draw}}-{{repID}}_two_sites-{type}.touchdown", type=["SUB","TIME"])
#     output:
#         sub=ROOT_dir+"/outputs/simulation/m0gtr/pred/stats/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}-two_sites-SUB.pickle",
#         time=ROOT_dir+"/outputs/simulation/m0gtr/pred/stats/{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}-two_sites-TIME.pickle",
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/m0gtr/pred/step_two_sites/",
#         pattern="'{geneID}-GTRG4-{omega}-{CpG}-{draw}-{repID}*'"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --output_sub={output.sub} \
#         --output_time={output.time} && rm $(echo {params.input_dir}{params.pattern})
#         """
