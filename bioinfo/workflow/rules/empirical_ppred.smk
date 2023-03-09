import os
import sys
import re
from pathlib import Path

ROOT_dir: str = Path(os.path.dirname(workflow.snakefile)).parent.__str__()
if ROOT_dir not in sys.path:
    sys.path.append(ROOT_dir)


def genererate_metadata(wildcards):
    s:str = '"{'
    geneID = wildcards.geneID
    s += "'geneID':'"+ geneID+"'"
    if "repID" in wildcards:
        repID = wildcards.repID
        s+= ",'repID':'"+repID+"'"
    s+='"}'
    return s



# rule fasta2phylip:
#     input:
#         script=ROOT_dir+"/scripts/data_prep.py",
#         data=ROOT_dir+"/data/{geneID}.fasta"
#     output: ROOT_dir+"/outputs/empirical/data/{geneID}.phylip"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --fasta={input.data} \
#         --phylip={output}
#         """

# rule generate_pb_mpi_cmd:
#     input:
#         script=ROOT_dir+"/scripts/generate_GTRG4_pb_mpi_cmd.py",
#         phylip=ROOT_dir+"/outputs/empirical/data/{geneID}.phylip",
#         tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
#     output:
#         docker=ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.sh",
#         cluster=ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}-cluster.sh"
#     params:
#         local=ROOT_dir,
#         docker="/data",
#         phylip="/outputs/empirical/data/{geneID}.phylip",
#         tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
#         cluster="/home/sll/workflow",
#         chainname="{geneID}-GTRG4-{repID}",
#         np="4",
#         output="{geneID}-GTRG4-{repID}-cluster.sh"
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --phylip={input.phylip} \
#         --tree={input.tree} \
#         --output={output.docker} \
#         --local={params.local} \
#         --docker={params.docker}
#         echo #####
#         echo "#qsub -q all.q@compute-0-1.local -q all.q@compute-0-2.local -q all.q@compute-0-3.local -q all.q@compute-0-4.local -q all.q@compute-0-5.local -q all.q@compute-0-6.local -q all.q@compute-0-7.local -q all.q@compute-0-0.local -cwd -pe mpi {params.np} {params.output}" >> {output.cluster}
#         echo "LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH" >> {output.cluster}
#         echo "CXX=`which g++`" >> {output.cluster}
#         echo "CC=`which gcc`" >> {output.cluster}
#         echo "export LD_LIBRARY_PATH" >> {output.cluster}
#         echo "export CC" >> {output.cluster}
#         echo "export CXX" >> {output.cluster}
#         echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} -gtr -ncat 1 -dgam 4 -s -x 10 200 {params.chainname} >> {output.cluster}
#         """



rule generate_M0GTR_pb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_M0GTRW_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/empirical/data/{geneID}.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:
        docker=ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.sh",
        cluster=ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}-cluster.sh"
    params:
        local=ROOT_dir,
        docker="/data",
        phylip="/outputs/empirical/data/{geneID}.phylip",
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
        echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} -mutsel -catfix uniform -freeomega -s -x 10 200 {params.chainname} >> {output.cluster}
        """


# rule generate_readpb_mpi_cmd:
#     input:
#         script=ROOT_dir+"/scripts/generate_readpb_mpi_ppred_cmd.py",
#         data=ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.chain"
#     output:
#         sh=ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}.sh",
#     params:
#         local=ROOT_dir,
#         docker="/data",
#         sampling="'100 1 200'",
#         model="GTRG4",
#         image="'ubuntu20.04/pbmpi:latest'",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input={input.data} \
#         --output={output.sh} \
#         --model={params.model} \
#         --image={params.image} \
#         --local={params.local} \
#         --docker={params.docker} \
#         --sampling={params.sampling} \
#         """


# rule run_NT_readpb_mpi_cmd:
#     input:
#         sh=ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}.sh",
#         phylip=ROOT_dir+"/outputs/empirical/data/{geneID}.phylip",
#         param=ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.param",
#     output:
#         touchdown=ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}_ppred.touchdown",
#         output_2=ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}_ppred.zip"
#     params:
#         sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/empirical/data/{geneID}.phylip")+"/'",
#         zip_distat_cmd="zip " +  ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}_ppred.zip " +  ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}_ppred*.ali",
#         mv_distat_cmd="mv "+ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}_ppred.zip "+ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/",
#         rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/empirical/pbmpi_gtrg4/{geneID}-GTRG4-{repID}_*.ali",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     threads: 2
#     shell:
#         """
#         {params.sed_cmd} {input.param} && bash {input.sh} && \
#         {params.zip_distat_cmd} && {params.rm_distat_s_cmd} && {params.mv_distat_cmd} && \
#         touch {output.touchdown}
#         """


rule generate_readpb_mpi_cmd_:
    input:
        script=ROOT_dir+"/scripts/generate_readpb_mpi_ppred_cmd.py",
        data=ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.chain"
    output:
        sh=ROOT_dir+"/outputs/empirical/ppred_test/readpb_m0gtr/{geneID}-M0GTR-{repID}.sh",
    params:
        local=ROOT_dir,
        docker="/data",
        sampling="'100 1 200'",
        model="M0GTR",
        image="'ubuntu20.04/pbmpi:latest'",
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


rule run_NT_readpb_mpi_cmd:
    input:
        sh=ROOT_dir+"/outputs/empirical/ppred_test/readpb_m0gtr/{geneID}-M0GTR-{repID}.sh",
        phylip=ROOT_dir+"/outputs/empirical/data/{geneID}.phylip",
        param=ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}.param",
    output:
        touchdown=ROOT_dir+"/outputs/empirical/ppred_test/readpb_m0gtr/{geneID}-M0GTR-{repID}_ppred.touchdown",
        output_2=ROOT_dir+"/outputs/empirical/ppred_test/readpb_m0gtr/{geneID}-M0GTR-{repID}_ppred.zip"
    params:
        sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/empirical/data/{geneID}.phylip")+"/'",
        zip_distat_cmd="zip " +  ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}_ppred.zip " +  ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}_ppred*.ali",
        mv_distat_cmd="mv "+ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}_ppred.zip "+ROOT_dir+"/outputs/empirical/ppred_test/readpb_m0gtr/",
        rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/empirical/pbmpi_m0gtr/{geneID}-M0GTR-{repID}_*.ali",
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
        zip=ROOT_dir+"/outputs/empirical/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}_ppred.zip"
    output:
        ROOT_dir+"/outputs/empirical/ppred_test/stats/{geneID}-GTRG4-{repID}_ppred.tsv"
    params:
        metadata=genererate_metadata,
    conda:
        ROOT_dir+"/envs/env.yml"
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.zip} --output={output} --metadata={params.metadata}
        """


rule compute_CpGfreq_:
    input:
        script=ROOT_dir+"/scripts/compute_CpGfreq.py",
        zip=ROOT_dir+"/outputs/empirical/ppred_test/readpb_m0gtr/{geneID}-M0GTR-{repID}_ppred.zip"
    output:
        ROOT_dir+"/outputs/empirical/ppred_test/stats/{geneID}-M0GTR-{repID}_ppred.tsv"
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
        ROOT_dir+"/outputs/empirical/ppred_test/stats/{geneID}-1_0-OBSERVED.tsv"
    params:
        metadata=genererate_metadata,
    conda:
        ROOT_dir+"/envs/env.yml"
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.fasta} --output={output} --metadata={params.metadata}
        """


# rule compute_stats:
#     input:
#         script=ROOT_dir+"/scripts/compute_CpG_stat_.py",
#         zipf=ROOT_dir+"/outputs/empirical/readpb_gtrg4/{geneID}-GTRG4-{repID}_suffdistatmap.zip",
#         prev_touchdown=ROOT_dir+"/outputs/empirical/readpb_gtrg4/{geneID}-GTRG4-{repID}_suffstatmap.touchdown",
#     output:
#         ROOT_dir+"/outputs/empirical/stats_gtrg4/{geneID}-GTRG4-{repID}_suffdistatmap_DINT.tsv",
#     params:
#         metadata=genererate_metadata
#     threads: 1
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} --input={input.zipf} --output={output} --metadata={params.metadata} && touch {output}
#         """


# rule compute_stats_:
#     input:
#         script=ROOT_dir+"/scripts/compute_CpG_stat_.py",
#         zipf=ROOT_dir+"/outputs/empirical/readpb_m0gtr/{geneID}-GTRG4-{repID}_suffdistatmap.zip",
#         prev_touchdown=ROOT_dir+"/outputs/empirical/readpb_gtrg4/{geneID}-GTRG4-{repID}_suffstatmap.touchdown",
#     output:
#         ROOT_dir+"/outputs/empirical/stats_gtrg4/{geneID}-GTRG4-{repID}_suffdistatmap_DINT.tsv",
#     params:
#         metadata=genererate_metadata
#     threads: 1
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} --input={input.zipf} --output={output} --metadata={params.metadata} && touch {output}
#         """
