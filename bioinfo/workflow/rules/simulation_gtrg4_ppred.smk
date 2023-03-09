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
    draw = wildcards.draw
    s += "'geneID':'"+ geneID+"','repID':'"+repID+"','draw':'"+draw+"'"+'"}'
    return s

rule fasta2phylip:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/data/{geneID}.fasta"
    output:
        ROOT_dir+"/outputs/simulation/gtrg4/data/{geneID}.phylip"
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} \
        --fasta={input.data} \
        --phylip={output}
        """

rule generate_GTRG4_pb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_GTR_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/simulation/gtrg4/data/{geneID}.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:
        docker=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.sh",
        cluster=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4/{geneID}-GTRG4-{repID}-cluster.sh"
    params:
        local=ROOT_dir,
        docker="/data",
        phylip="/outputs/simulation/gtrg4/data/{geneID}.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        cluster="/home/sll/mappings_brief",
        chainname="{geneID}-GTRG4-{repID}",
        np="4",
        output="{geneID}-GTRG4-{repID}-cluster.sh"
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
        echo "#qsub -cwd -pe mpi {params.np} {params.output}" >> {output.cluster}
        echo "LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH" >> {output.cluster}
        echo "CXX=`which g++`" >> {output.cluster}
        echo "CC=`which gcc`" >> {output.cluster}
        echo "export LD_LIBRARY_PATH" >> {output.cluster}
        echo "export CC" >> {output.cluster}
        echo "export CXX" >> {output.cluster}
        echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} -gtrg4 -ncat 1 -dgam 1 -s -x 1 2000 {params.chainname} >> {output.cluster}
        """

rule generate_gtrg4_values:
    input:
        script=ROOT_dir+"/scripts/generate_GTRG4_values.py",
        mcmc=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.chain",
    output:
        values=ROOT_dir+"/outputs/simulation/gtrg4/values/{geneID}-GTRG4-{repID}-{draw}.values",
        tree=ROOT_dir+"/outputs/simulation/gtrg4/values/{geneID}-GTRG4-{repID}-{draw}.tre"
    params:
        burnin=1000,
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

rule generate_gtrg4_ppred_values:
    input:
        script=ROOT_dir+"/scripts/generate_GTRG4_values.py",
        mcmc=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4/{geneID}-GTRG4-{repID}.chain",
    output:
        values=ROOT_dir+"/outputs/simulation/gtrg4/ppred/values/{geneID}-GTRG4-{repID}-{draw}.values",
        tree=ROOT_dir+"/outputs/simulation/gtrg4/ppred/values/{geneID}-GTRG4-{repID}-{draw}.tre"
    params:
        burnin=1000,
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


def get_Nsite(wildcards):
    geneID = wildcards.geneID
    return int(config["Nsites"][geneID])

rule generate_simu_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_GTR_simu_cmd.py",
        values=ROOT_dir+"/outputs/simulation/gtrg4/values/{geneID}-GTRG4-{repID}-{draw}.values",
        tree=ROOT_dir+"/outputs/simulation/gtrg4/values/{geneID}-GTRG4-{repID}-{draw}.tre"
    output:
        sh=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.sh"
    params:
        simu=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}",
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

rule generate_ppred_simu_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_GTR_simu_cmd.py",
        values=ROOT_dir+"/outputs/simulation/gtrg4/ppred/values/{geneID}-GTRG4-{repID}-{draw}.values",
        tree=ROOT_dir+"/outputs/simulation/gtrg4/ppred/values/{geneID}-GTRG4-{repID}-{draw}.tre"
    output:
        sh=ROOT_dir+"/outputs/simulation/gtrg4/ppred/simu/{geneID}-GTRG4-{repID}-{draw}.sh"
    params:
        simu=ROOT_dir+"/outputs/simulation/gtrg4/ppred/simu/{geneID}-GTRG4-{repID}-{draw}",
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
        ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.sh"
    output:
        simu=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.fa",
    threads: 1
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule run_ppred_simu_cmd:
    input:
        ROOT_dir+"/outputs/simulation/gtrg4/ppred/simu/{geneID}-GTRG4-{repID}-{draw}.sh"
    output:
        simu=ROOT_dir+"/outputs/simulation/gtrg4/ppred/simu/{geneID}-GTRG4-{repID}-{draw}.fa",
    threads: 1
    shell:
        """
        echo bash {input}
        bash {input}
        """

rule fasta2phylip_simu:
    input:
        script=ROOT_dir+"/scripts/data_prep.py",
        data=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.fa"
    output:
        ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip"
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
        fasta=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.fa"
    output:
        ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/stats/{geneID}-GTRG4-{repID}-{draw}-OBSERVED.tsv"

    conda:
        ROOT_dir+"/envs/env.yml"

    params:
        metadata=genererate_metadata
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.fasta} --output={output} --metadata={params.metadata}
        """

rule generate_GTRG4__pb_mpi_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_GTRG4_pb_mpi_cmd.py",
        phylip=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip",
        tree=ROOT_dir+"/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    output:
        docker=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.sh",
        cluster=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}-cluster.sh"
    params:
        local=ROOT_dir,
        docker="/data",
        phylip="/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip",
        tree="/data/Mammalia-39sp-CAT_unrooted_withoutsupport.tre",
        cluster="/home/sll/mappings_brief",
        chainname="{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}",
        np="4",
        output="{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}-cluster.sh"
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
        echo "#qsub -cwd -pe mpi {params.np} {params.output}" >> {output.cluster}
        echo "LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH" >> {output.cluster}
        echo "CXX=`which g++`" >> {output.cluster}
        echo "CC=`which gcc`" >> {output.cluster}
        echo "export LD_LIBRARY_PATH" >> {output.cluster}
        echo "export CC" >> {output.cluster}
        echo "export CXX" >> {output.cluster}
        echo mpirun -np {params.np} /home/sll/pbmpi/data/pb_mpi -d {params.cluster}{params.phylip} -T {params.cluster}{params.tree} -gtr -ncat 1 -dgam 4 -s -x 1 2000 {params.chainname} >> {output.cluster}
        """

rule generate_readpb_mpi_ppred_cmd:
    input:
        script=ROOT_dir+"/scripts/generate_readpb_mpi_ppred_cmd.py",
        data=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.chain"
    output:
        sh=ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.sh",
    params:
        local=ROOT_dir,
        docker="/data",
        sampling="'1000 20 2000'",
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
        sh=ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.sh",
        phylip=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip",
        param=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.param",
    output:
        touchdown=ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.touchdown",
        output_2=ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.zip"
    params:
        sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip")+"/'",
        zip_distat_cmd="zip " +  ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.zip " +  ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.ali",
        mv_distat_cmd="mv "+ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.zip "+ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/",
        rm_distat_s_cmd="rm "+ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.ali",

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
        zip=ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.zip"
    output:
        ROOT_dir+"/outputs/simulation/gtrg4/ppred_test/stats/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_ppred.tsv"
    conda:
        ROOT_dir+"/envs/env.yml"
    params:
        metadata=genererate_metadata
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.zip} --output={output} --metadata={params.metadata}
        """



rule compute_stats_DINT:
    input:
        script=ROOT_dir+"/scripts/compute_CpG_stat_.py",
        zipf=ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_suffdistatmap.zip"
    output:
        ROOT_dir+"/outputs/simulation/gtrg4/stats_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_suffdistatmap_DINT.tsv",
    params:
        metadata=genererate_metadata
    threads: 1
    shell:
        """
        python3 {input.script} --input={input.zipf} --output={output} --metadata={params.metadata} && touch {output}
        """


rule aggregated_stats_DINT:
    input:
        script=ROOT_dir+"/scripts/aggregate_CpG_stats.py",
    output:
         ROOT_dir+"/outputs/simulation/gtrg4/stats_gtrg4/aggregated_stats_DINT.tsv"
    params:
        input_dir=ROOT_dir+"/outputs/simulation/gtrg4/stats_gtrg4/",
        pattern="*_suffdistatmap_DINT.tsv",
    conda:
        ROOT_dir+"/envs/env.yml"
    shell:
        """
        python3 {input.script} --input_dir={params.input_dir} --pattern={params.pattern} --output={output}
        """





# rule run_readpb_mpi_cmd:
#     input:
#         sh=ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.sh",
#         phylip=ROOT_dir+"/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip",
#         param=ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}.param",
#     output:
#         ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_map.touchdown"
#     params:
#         sed_cmd="sed -i '7s/.*/"+re.sub("/","\/","/data/outputs/simulation/gtrg4/simu/{geneID}-GTRG4-{repID}-{draw}.phylip")+"/'",
#         rm_cmd_nodestates=get_nodestates_to_delete,
#         mv_cmd="mv "+ROOT_dir+"/outputs/simulation/gtrg4/pbmpi_gtrg4_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.map "+ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/",
#     threads: 2
#     shell:
#         """
#         {params.sed_cmd} {input.param} && bash {input.sh} && {params.rm_cmd_nodestates} && {params.mv_cmd} && touch {output}
#         """

# rule serialize_mapping:
#     input:
#         script=ROOT_dir+"/scripts/serialize_mapping_.py",
#     output:
#         ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.map'",
#         output_dir=ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/",
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
# #&& rm -f $(echo {params.input_dir}{params.pattern})
# rule generate_two_sites_mappings:
#     input:
#         script=ROOT_dir+"/scripts/generate_two_sites_mapping_.py",
#         prevtouchdown=ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites_/",
#         list_of_paired_sites=get_adjsites
#     output:
#         ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites.touchdown"
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
#         prevtouchdown=ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites_/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites/",
#     output:
#         sub=ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites-SUB.touchdown",
#         time=ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites-TIME.touchdown",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --output_dir={params.output_dir}  && touch {output.sub} && touch {output.time}
#         """
# # && rm -f $(echo {params.input_dir}{params.pattern})

# rule aggregate_two_sites_results:
#     input:
#         script=ROOT_dir+"/scripts/aggregate_results_over_sites.py",
#         prevtouchdown=expand(ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites/{{geneID}}-GTRG4-{{repID}}-{{draw}}-GTRG4-{{repID_}}_two_sites-{type}.touchdown", type=["SUB","TIME"])
#     output:
#         sub=ROOT_dir+"/outputs/simulation/gtrg4/stats/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}-two_sites-SUB.pickle",
#         time=ROOT_dir+"/outputs/simulation/gtrg4/stats/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}-two_sites-TIME.pickle",
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/step_two_sites/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}*'"
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
# #&& rm $(echo {params.input_dir}{params.pattern})

# rule serialize_mapping_pred:
#     input:
#         script=ROOT_dir+"/scripts/serialize_mapping_.py",
#     output:
#         ROOT_dir+"/outputs/simulation/gtrg4/pred/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.map'",
#         output_dir=ROOT_dir+"/outputs/simulation/gtrg4/pred/readpb_gtrg4/",
#         metadata=get_params_pred,
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
#         #&& rm -f $(echo {params.input_dir}{params.pattern})

# rule generate_two_sites_mappings_pred:
#     input:
#         script=ROOT_dir+"/scripts/generate_two_sites_mapping_.py",
#         prevtouchdown=ROOT_dir+"/outputs/simulation/gtrg4/pred/readpb_gtrg4/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_pickle.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/pred/readpb_gtrg4/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites_/",
#         list_of_paired_sites=get_adjsites
#     output:
#         ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites.touchdown"
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
#         prevtouchdown=ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites_/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites.touchdown"
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites_/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_*.pickle'",
#         output_dir=ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites/",
#     output:
#         sub=ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites-SUB.touchdown",
#         time=ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}_two_sites-TIME.touchdown",
#     conda:
#         ROOT_dir+"/envs/env.yml"
#     shell:
#         """
#         python3 {input.script} \
#         --input_dir={params.input_dir} \
#         --pattern={params.pattern} \
#         --output_dir={params.output_dir}  && touch {output.sub} && touch {output.time}
#         """
# #&& rm -f $(echo {params.input_dir}{params.pattern})

# rule aggregate_two_sites_results_pred:
#     input:
#         script=ROOT_dir+"/scripts/aggregate_results_over_sites.py",
#         prevtouchdown=expand(ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites/{{geneID}}-GTRG4-{{repID}}-{{draw}}-GTRG4-{{repID_}}_two_sites-{type}.touchdown", type=["SUB","TIME"])
#     output:
#         sub=ROOT_dir+"/outputs/simulation/gtrg4/pred/stats/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}-two_sites-SUB.pickle",
#         time=ROOT_dir+"/outputs/simulation/gtrg4/pred/stats/{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}-two_sites-TIME.pickle",
#     params:
#         input_dir=ROOT_dir+"/outputs/simulation/gtrg4/pred/step_two_sites/",
#         pattern="'{geneID}-GTRG4-{repID}-{draw}-GTRG4-{repID_}*'"
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
# #&& rm $(echo {params.input_dir}{params.pattern})
