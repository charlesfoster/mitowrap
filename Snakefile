#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 15:43:06 2021

@author: cfos
"""
###########
# Libraries
###########

import os
import glob
import re
from datetime import date
import pandas as pd
import yaml
import glob
from Bio import SeqIO
import subprocess

#################
# Custom functions
#################

def find_input_data(READSDIR, data_type, SUFFIX):
    reads = glob.glob(os.path.join(READSDIR,'*'+SUFFIX))
    main_exclusion = ['Undet','NEG_', 'NC_']
    neg_inclusion = ['NEG_', 'NC_']
    if data_type == "main":
        reads = [i for i in reads if not any(b in i for b in main_exclusion)]
    elif data_type == "neg":
        reads = [i for i in reads if any(b in i for b in neg_inclusion)]
    result = []
    for read in reads:
        if SUFFIX == '_L001_R1_001.fastq.gz':
            sample = re.sub("_S\d+_L001.*", repl = "", string=os.path.basename(read))
            barcode = re.search("S\d+_L001.*", string=os.path.basename(read)).group()
            barcode = re.sub("_.*",repl="",string=barcode)
            res = {'sample':sample,'barcode':barcode, 'fq1':read,'fq2':read.replace("_R1","_R2")}
        else:
            sample = re.sub(SUFFIX, repl = "", string=os.path.basename(read))
            barcode = "N/A"
            res = {'sample':sample,'barcode':barcode, 'fq1':read,'fq2':read.replace("_R1","_R2")}    
        result.append(res)
    result =  pd.DataFrame(result, dtype=str).set_index(["sample"], drop=False)
    return result


def get_trim_names(wildcards):
    """
    This function:
      1. Returns the correct input and output trimmed file names for fastp.
    """
    inFile = INPUT_TABLE.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    return "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " + os.path.join(RESULT_DIR, wildcards.sample, "fastp", wildcards.sample + "_trimmed_R1.fq.gz") + " --out2 " + os.path.join(RESULT_DIR, wildcards.sample, "fastp", wildcards.sample + "_trimmed_R2.fq.gz")

def get_fastq(wildcards):
    """This function returns the forward and reverse fastq files for samples"""
    return INPUT_TABLE.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()


###############
# Configuration
###############
configfile: "config.yaml"

onstart:
    header = """
----------------------------------------------
\033[95m
           _ _                                
 _ __ ___ (_) |_ _____      ___ __ __ _ _ __  
| '_ ` _ \| | __/ _ \ \ /\ / / '__/ _` | '_ \ 
| | | | | | | || (_) \ V  V /| | | (_| | |_) |
|_| |_| |_|_|\__\___/ \_/\_/ |_|  \__,_| .__/ 
                                       |_|        
\033[0m
----------------------------------------------

"""
    print(header)
    print("~~~ CONFIGURATION ~~~")
    for item in config:
        print(item+": "+str(config[item]))
    print("\n")
main_dir = os.path.abspath(os.path.dirname("config.yaml"))
REFERENCE = config['reference']
READSDIR = config["reads_dir"]
RESULT_DIR = config["outdir"]

if not os.path.exists(RESULT_DIR):
    os.makedirs(RESULT_DIR)
SUFFIX = config['suffix']
SUFFIX = SUFFIX.replace('R2','R1')
KEEP_READS = config["keep_reads"]

MAIN_TABLE = find_input_data(READSDIR, "main", SUFFIX)
MAIN_TABLE.to_csv( os.path.join(RESULT_DIR,"sample_sheet.csv"), index=False, header=True)
MAIN_SAMPLES = MAIN_TABLE.index.get_level_values('sample').unique().tolist()
INPUT_TABLE = find_input_data(READSDIR, "all", SUFFIX)

#if ISOLATES != False:
#    MAIN_SAMPLES = [i for i in MAIN_SAMPLES if any(b in i for b in ISOLATES)]

config['analysis_samples'] = '; '.join(MAIN_SAMPLES)

TODAY=date.today().strftime("%Y-%m-%d")

if config['using_conda']:
    os.system('echo "True" > {}'.format(os.path.join(main_dir, ".using_conda")))
    using_conda = "True"
else:
    using_conda = "False"
    os.system('echo "False" > {}'.format(os.path.join(main_dir, ".using_conda")))

if os.path.exists(os.path.join(main_dir, ".ete_data_added")):
    os.remove(os.path.join(main_dir, ".ete_data_added"))

################
# Optional removal of trimmed reads (to save space)
################
if KEEP_READS == False:
    onsuccess:
        print("Removing trimmed reads (input reads remain untouched)\n")
        dead_reads = glob.glob(RESULT_DIR + '/**/*trimmed*.gz', recursive=True)
        [os.remove(x) for x in dead_reads]


################
# Desired outputs
################

rule final_qc:
    input:
        sample_reports = expand(os.path.join(RESULT_DIR, "{sample}.qc_results.csv"), sample = MAIN_SAMPLES)
    message:
        "Pipeline complete!"
    run:
        qc_files = glob.glob(RESULT_DIR+'/*qc_results.csv')
        combined_csv = pd.concat([pd.read_csv(f) for f in qc_files ])
        combined_csv.to_csv( os.path.join(RESULT_DIR,TODAY+"_QC.csv"), index=False, header=True)
        with open(os.path.join(RESULT_DIR,"config.yaml"), 'w') as outfile:
            yaml.dump(config, outfile, default_flow_style=False)


########################
# Rules 
########################

rule fastp:
    input:
        get_fastq,
    output:
        fq1  = os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R1.fq.gz"),
        fq2  = os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R2.fq.gz"),
        html = os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_fastp.html"),
        json = os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_fastp.json")
    message:"trimming {wildcards.sample} reads"
    threads: 4
    log:
        os.path.join(RESULT_DIR,"{sample}/fastp/{sample}.log.txt")
    params:
        in_and_out_files =  get_trim_names,
        sampleName = "{sample}",
        qualified_quality_phred = 20,
        cut_mean_quality = 20,
        unqualified_percent_limit = 10,
        length_required = 50
    container:
        "docker://quay.io/biocontainers/fastp:0.23.2--hb7a2d85_2"
    conda:
        "./envs/fastp.yaml"
    resources: cpus=4
    shell:
        """
        fastp --thread {threads} \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_mean_quality 20 \
        --qualified_quality_phred {params.qualified_quality_phred} \
        --unqualified_percent_limit {params.unqualified_percent_limit} \
        --correction \
        --length_required 50 \
        --html {output.html} \
        --json {output.json} \
        {params.in_and_out_files} \
        2>{log}
        """

rule bwa_index:
    input:
        reference = REFERENCE,
    output:
        index = REFERENCE+".bwt",
    threads: 1
    resources:
        cpus=1,
    container:
        "docker://biocontainers/bwa:v0.7.17_cv1"
    conda:
        "./envs/bwa.yaml"
    shell:
        """
        bwa index {input.reference}
        """

rule bwa_map:
    input:
        r1=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R1.fq.gz"),
        r2=os.path.join(RESULT_DIR, "{sample}/fastp/{sample}_trimmed_R2.fq.gz"),
        index = REFERENCE+".bwt",
    output:
        sam=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.sam")),
    message:
        "mapping {wildcards.sample} reads to reference"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/{sample}.bwa.log"),
    params:
        reference=REFERENCE,
    resources:
        cpus=4,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    container:
        "docker://biocontainers/bwa:v0.7.17_cv1"
    conda:
        "./envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {params.reference} \
        {input.r1} {input.r2} > {output.sam} 2>{log}
        """


rule samtools_fastq:
    input:
        sam = os.path.join(RESULT_DIR, "{sample}/{sample}.sam"),
    output:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.mapped.bam"),
        r1_mapped=os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R1.fq.gz"),
        r2_mapped=os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R2.fq.gz"),
        singletons_mapped=os.path.join(RESULT_DIR, "{sample}/{sample}_mito_singletons.fq.gz"),
        r1_unmapped=os.path.join(RESULT_DIR, "{sample}/{sample}_unmapped_R1.fq.gz"),
        r2_unmapped=os.path.join(RESULT_DIR, "{sample}/{sample}_unmapped_R2.fq.gz"),
        singletons_unmapped=os.path.join(RESULT_DIR, "{sample}/{sample}_unmapped_singletons.fq.gz"),
    message:
        "extracting reads for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/{sample}.samtools.log"),
    params:
        reference=REFERENCE,
    resources:
        cpus=4,
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*",
    container:
        "docker://staphb/samtools:1.15"
    conda:
        "./envs/samtools.yaml"
    shell:
        """
        samtools view -@ {threads} -u {input.sam} 2> {log} | \
        samtools sort -n -@ {threads} -o {output.bam} - 2>> {log} 
        samtools fastq -F 4 -c 6 -@ {threads} -1 {output.r1_mapped} -2 {output.r2_mapped} -s {output.singletons_mapped} -N {output.bam} 2>> {log}  
        samtools fastq -f 4 -c 6 -@ {threads} -1 {output.r1_unmapped} -2 {output.r2_unmapped} -s {output.singletons_unmapped} -N {output.bam} 2>> {log} 
        """


rule add_animal_mt_db:
    output:
        db_added = os.path.join(main_dir, ".animal_db_added"),
    threads: 1
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*"
    container:
        "docker://quay.io/biocontainers/getorganelle:1.7.6.1--pyh5e36f6f_0",
    log:
        ".animal_db.log"
    conda:
        "./envs/get_organelle.yaml",
    shell:
        """
        get_organelle_config.py -a animal_mt >> {log} 2>&1
        touch {output.db_added}
        """


rule get_organelle_assembly:
    input:
        db_added = os.path.join(main_dir, ".animal_db_added"),
        r1 = os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R1.fq.gz"),
        r2 = os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R2.fq.gz")
    output:
        ckp = os.path.join(RESULT_DIR, "{sample}/assembly/{sample}.getOrgComplete.txt"),
        complete_file = os.path.join(RESULT_DIR, "{sample}/assembly/get_org.log.txt")
    params: 
        sample = "{sample}",
        outdir = os.path.join(RESULT_DIR, "{sample}/assembly")
    params: 
    message:"assembling {wildcards.sample} using getOrganelle"
    threads: 8
    log:
        os.path.join(RESULT_DIR,"{sample}/{sample}.getOrg.log")
    resources:
        mem_mb=16000,
        cpus=8
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*"
    container:
        "docker://quay.io/biocontainers/getorganelle:1.7.6.1--pyh5e36f6f_0",
    conda:
        "./envs/get_organelle.yaml",
    shell:
        """
        echo "Using getOrganelle" >> {log} 2>&1
        get_organelle_from_reads.py -1 {input.r1} -2 {input.r2} \
        -R 10 -k 21,45,65,85,105 -F animal_mt -o {params.outdir} -t {threads} --overwrite >> {log} 2>&1
        touch {output.ckp}
        """


rule get_etetoolkit_data:
    input:
        using_conda = os.path.join(main_dir, ".using_conda"),
    output:
        ete_data = os.path.join(main_dir, ".ete_data_added"),
    container:
        "docker://guanliangmeng/mitoz:3.4"
    message:"Checking for ete3 NCBI taxonomy"
    conda:
        "./envs/mitoz.yaml"
    shell:
        "python3 ./scripts/get_etetoolkit_data.py {input.using_conda} {output.ete_data}"

rule mitoz_assembly:
    input:
        ete_data = os.path.join(main_dir, ".ete_data_added"),
        r1 = os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R1.fq.gz"),
        r2 = os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R2.fq.gz"),

    output:
        ckp = os.path.join(RESULT_DIR, "{sample}/assembly/{sample}.mitozComplete.txt"),
        complete_file = os.path.join(RESULT_DIR, "{sample}/mitoz/{sample}.result/{sample}.{sample}.megahit.mitogenome.fa.result/summary.txt"),
        mitoz_assembly = os.path.join(RESULT_DIR, "{sample}/mitoz/{sample}.result/{sample}.megahit.result/{sample}.megahit.mitogenome.fa"),
    params: 
        sample = "{sample}",
        clade = config["clade"],
        species_name = config["species_name"],
        outdir = os.path.join(RESULT_DIR, "{sample}/assembly"),
        workdir = os.path.join(RESULT_DIR, "{sample}","mitoz"),
        touchdir1 = os.path.join(RESULT_DIR, "{sample}/mitoz/{sample}.result/{sample}.{sample}.megahit.result/"),
        touchdir2 = os.path.join(RESULT_DIR, "{sample}/mitoz/{sample}.result/{sample}.{sample}.megahit.mitogenome.fa.result/"),
    message:"assembling {wildcards.sample} using MitoZ"
    threads: 8
    log:
        os.path.join(RESULT_DIR,"{sample}/{sample}.mitoz.log")
    resources:
        mem_mb=30000,
        cpus=8
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*"
    container:
        "docker://guanliangmeng/mitoz:3.4"
    conda:
        "./envs/mitoz.yaml"
    shell:
        """
        mitoz all \
        --workdir {params.workdir} \
        --skip_filter \
        --outprefix {wildcards.sample} \
        --thread_number 8 \
        --clade {params.clade} \
        --genetic_code auto \
        --species_name {params.species_name} \
        --fq1 {input.r1} \
        --fq2 {input.r2} \
        --fastq_read_length 150 \
        --data_size_for_mt_assembly 0 \
        --assembler megahit \
        --kmers_megahit 39 59 79 99 119 141 \
        --memory 50 \
        --requiring_taxa {params.clade} >> {log} 2>&1 || mkdir -p {params.touchdir1} && mkdir -p {params.touchdir2} && touch {output.complete_file} && touch {output.mitoz_assembly}

        touch {output.ckp}
        """

rule annotate_getOrg_assembly:
    input:
        r1 = os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R1.fq.gz"),
        r2 = os.path.join(RESULT_DIR, "{sample}/{sample}_mito_R2.fq.gz"),
        ckp = os.path.join(RESULT_DIR, "{sample}/assembly/{sample}.getOrgComplete.txt"),
    output:
        ckp = os.path.join(RESULT_DIR, "{sample}/assembly/{sample}.annotateGetOrgComplete.txt"),
        summary = os.path.join(RESULT_DIR, "{sample}","getOrg_annotation","{sample}.{sample}.assembly.fa.result","summary.txt"),
    params: 
        clade = config["clade"],
        species_name = config["species_name"],
        assembly_dir = os.path.join(RESULT_DIR, "{sample}","assembly"),
        annotation_dir = os.path.join(RESULT_DIR, "{sample}","getOrg_annotation"),
        new_fasta = os.path.join(RESULT_DIR, "{sample}","getOrg_annotation","{sample}.assembly.fa")
    message:"annotating getOrganelle assembly for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR,"{sample}/{sample}.mitozAnnotateGetOrg.log")
    resources:
        mem_mb=16000,
        cpus=4
    wildcard_constraints:
        sample="(?!NC)(?!NEG).*"
    container:
        "docker://guanliangmeng/mitoz:3.4"
    conda:
        "./envs/mitoz.yaml"
    shell:
        """
        FASTA=$(find {params.assembly_dir} -maxdepth 1 -name "*.fasta")
        mkdir -p {params.annotation_dir}
        if [[ -z $FASTA ]]; then
            printf ">{wildcards.sample}\n" > {params.new_fasta}
            ASSDIR=$(dirname {params.assembly_dir})
            ANNDIR=$(dirname {params.assembly_dir})
            mkdir -p $ASSDIR
            mkdir -p $ANNDIR
            touch {output.ckp}
            touch {output.summary}
        else
            cp $FASTA {params.new_fasta}
            SAMPLE={wildcards.sample}
            sed -i "s|>.*|>$SAMPLE|g" {params.new_fasta}
            mitoz annotate \
            --workdir {params.annotation_dir} \
            --thread_number 4  \
            --outprefix {wildcards.sample} \
            --fastafiles {params.new_fasta} \
            --fq1 {input.r1} \
            --fq2 {input.r1} \
            --species_name {params.species_name} \
            --genetic_code auto \
            --clade {params.clade} >> {log} 2>&1
            touch {output.ckp}
        fi
        """

rule sample_qc:
    input:
        getOrg_summary = os.path.join(RESULT_DIR, "{sample}","getOrg_annotation","{sample}.{sample}.assembly.fa.result","summary.txt"),
        mitoz_summary = os.path.join(RESULT_DIR, "{sample}/mitoz/{sample}.result/{sample}.{sample}.megahit.mitogenome.fa.result/summary.txt"),
        mitoz_assembly = os.path.join(RESULT_DIR, "{sample}/mitoz/{sample}.result/{sample}.megahit.result/{sample}.megahit.mitogenome.fa"),
    params:
        joint_summary = os.path.join(RESULT_DIR, "{sample}/assembly_summaries.csv"),
        sample = "{sample}",
    output:
        report = temp(os.path.join(RESULT_DIR, "{sample}.qc_results.csv")),
    log:
        log = os.path.join(RESULT_DIR,"{sample}/{sample}.parse_summaries.log"),
    run:
        # parse the 'summary' files
        f1 = input.getOrg_summary
        f2 = input.mitoz_summary
        sample = params.sample
        header = 'sample_name,assembly_program,circular,total_PCG,total_tRNA,total_rRNA,total_genes'

        with open(params.joint_summary, 'w') as o:
            print(header,file=o)

        if os.path.getsize(f1) == 0:
            with open(params.joint_summary, 'a') as o:
                program = "getOrganelle"
                out = f"{sample},{program},no,0,0,0,0"
                print(out,file=o)
        else:
            with open(f1,'r') as infile:
                lines = infile.readlines()
            
            count = 0
            for line in lines:
                count += 1
                if count == 2:
                    circular = line.split()[2]
                if 'Protein coding genes totally found' in line:
                    total_pcg = line.split()[-1]
                if 'tRNA genes totally found' in line:
                    total_trna = line.split()[-1]
                if 'rRNA genes totally found' in line:
                    total_rrna = line.split()[-1]
                if 'Genes totally found' in line:
                    total_genes = line.split()[-1]
            with open(params.joint_summary, 'a') as o:
                out = f"{sample},getOrganelle,{circular},{total_pcg},{total_trna},{total_rrna},{total_genes}"
                print(out,file=o)

        if os.path.getsize(f2) == 0:
            with open(params.joint_summary, 'a') as o:
                program = "mitoz"
                out = f"{sample},{program},no,0,0,0,0"
                print(out,file=o)
        else:
            with open(f2,'r') as infile:
                lines = infile.readlines()
            
            count = 0
            for line in lines:
                count += 1
                if count == 2:
                    circular = line.split()[2]
                if 'Protein coding genes totally found' in line:
                    total_pcg = line.split()[-1]
                if 'tRNA genes totally found' in line:
                    total_trna = line.split()[-1]
                if 'rRNA genes totally found' in line:
                    total_rrna = line.split()[-1]
                if 'Genes totally found' in line:
                    total_genes = line.split()[-1]
            with open(params.joint_summary, 'a') as o:
                out = f"{sample},mitoz,{circular},{total_pcg},{total_trna},{total_rrna},{total_genes}"
                print(out,file=o)

        # parse the fasta files
        fastas = [os.path.join(RESULT_DIR, params.sample, "assembly", x) for x in os.listdir(os.path.join(RESULT_DIR, params.sample, 'assembly')) if x.endswith(".fasta")]
        fa_list = []
        if len(fastas) > 0:
            for fasta in fastas:
                fa_name = os.path.basename(fasta)
                for record in SeqIO.parse(fasta, "fasta"):
                    info = {'sample_name':params.sample,
                            'assembly_program':'getOrganelle', 
                            'filename':fa_name,
                            'contig':record.id,
                            'length':len(record.seq)}
                    fa_list.append(info)
        else:
            fa_list.append({'sample_name':params.sample,
            'assembly_program':'getOrganelle', 
            'filename':'Failed',
            'contig':None,
            'length':0})

        if os.path.getsize(input.mitoz_assembly) != 0:
            fa_name = os.path.basename(input.mitoz_assembly)
            for record in SeqIO.parse(input.mitoz_assembly, "fasta"):
                info = {'sample_name':params.sample,
                        'assembly_program':'mitoz', 
                        'filename':fa_name,
                        'contig':record.id,
                        'length':len(record.seq)}
                fa_list.append(info)
        else:
            fa_list.append({'sample_name':params.sample,
            'assembly_program':'mitoz', 
            'filename':'Failed',
            'contig':None,
            'length':0})

        # join everything
        result1 = pd.DataFrame.from_dict(fa_list)
        result2 = pd.read_csv(params.joint_summary)
        result = result1.merge(result2, on=['sample_name','assembly_program'], how='left')
        result.sort_values(by=['length'], inplace=True,ascending=False)
        result.to_csv(output.report, header=True,index=False)
