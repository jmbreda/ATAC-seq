import os
configfile: "config/male_female_liver_atac_seq.yaml"
#configfile: "config/test_atac_seq.yaml"
READS = ['1', '2']

rule all:
    input:
        expand(os.path.join(config['peak_dir'],'{sample}_accessible_regions.gappedPeak'),
            sample=config['SAMPLE']),
        expand(os.path.join(config['mapping_dir'],'{sample}_coverage.bw'),
            sample=config['SAMPLE'])
        #expand(os.path.join(config['peak_dir'],'{sample}.filteredPeaks.gappedPeak'),
        #    sample=config['SAMPLE']),
        #expand(os.path.join(config['peak_dir'],'{sample}.filteredSummits.bed'),
        #    sample=config['SAMPLE'])

rule fastqc:
    input:
        "resources/reads/{sample}_{read}.fastq.gz"
    output:
        os.path.join(config['fastqc_dir'],"{sample}_{read}_fastqc.html")
    log:
        "logs/fastqc/{sample}_{read}.log"
    threads: 1
    params:
        outdir=config['fastqc_dir']
    #conda:
    #	"../envs/fastqc.yaml"
    shell:
        """
        #!/bin/bash
        fastqc {input} --outdir={params.outdir}
        """

rule trimming:
    input:
        ["resources/reads/{{sample}}_{read}.fastq.gz".format(read=read) for read in READS]
    output:
        forward_paired = os.path.join(config['trimmomatic_dir'],"{sample}_forward_paired.fq.gz"),
        forward_unpaired = os.path.join(config['trimmomatic_dir'],"{sample}_forward_unpaired.fq.gz"),
        reverse_paired = os.path.join(config['trimmomatic_dir'],"{sample}_reverse_paired.fq.gz"),
        reverse_unpaired = os.path.join(config['trimmomatic_dir'],"{sample}_reverse_unpaired.fq.gz")
    #message: "Trimming Illumina adapters from {input.forward_read} and {input.reverse_read}"
    log:
        "logs/trimmomatic/{sample}.log"
    wildcard_constraints:
        sample='|'.join(config['SAMPLE'])
    shell:
        """
        trimmomatic PE {input} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
        ILLUMINACLIP:{config[adapter]} LEADING:{config[leading]} TRAILING:{config[trailing]} MINLEN:{config[minlen]}
        """

#rule bwa-mem2_index:
#	input:
#		config['genome']
#	shell:
#		"bwa-mem2 index {input}"

rule map:
    input:
        forward_paired = os.path.join(config['trimmomatic_dir'],"{sample}_forward_paired.fq.gz"),
        reverse_paired = os.path.join(config['trimmomatic_dir'],"{sample}_reverse_paired.fq.gz")
    params:
        index = config['genome']
    output:
        os.path.join(config['mapping_dir'],'{sample}.sam')
    threads: 12
    shell:
        "bwa-mem2 mem -t {threads} {params.index} {input.forward_paired} {input.reverse_paired} > {output}"

rule samtools_sort:
    input:
        os.path.join(config['mapping_dir'],'{sample}.sam')
    output:
        os.path.join(config['mapping_dir'],'{sample}_sorted.bam')
    threads: 4
    shell:
        "samtools sort {input} -o {output} -@ {threads}"

rule samtools_index:
    input:
        os.path.join(config['mapping_dir'],'{sample}_sorted.bam')
    output:
        os.path.join(config['mapping_dir'],'{sample}_sorted.bam.bai')
    threads: 4
    shell:
        "samtools index -@ {threads} {input} {output}"

rule samtools_genome_info:
    input:
        os.path.join(config['mapping_dir'],'{sample}_sorted.bam')
    output:
        os.path.join(config['mapping_dir'],'{sample}_genome.info')
    shell:
        "scripts/samtools_genome_info.sh {input} {output}"

rule coverage:
    input:
        os.path.join(config['mapping_dir'],'{sample}_sorted.bam')
    output:
        os.path.join(config['mapping_dir'],'{sample}_coverage.bw')
    shell:
        """
        module load gcc/11.3.0
        module load py-deeptools
        bamCoverage -b {input} -o {output}
        """

rule peak_calling:
    resources:
        nodes=1,
        task=1,
        cpu_per_task=1,
        mem="200G",
        time="12:00:00"
    input:
        bam = os.path.join(config['mapping_dir'],'{sample}_sorted.bam'),
        index = os.path.join(config['mapping_dir'],'{sample}_sorted.bam.bai'),
        genome = os.path.join(config['mapping_dir'],'{sample}_genome.info')
    output:
        #peak=os.path.join(config['peak_dir'],'{sample}_peaks.gappedPeak'),
        peak=os.path.join(config['peak_dir'],'{sample}_accessible_regions.gappedPeak')
        #summit=os.path.join(config['peak_dir'],'{sample}_summits.bed')
    #conda: "/home/jbreda/M-F_Liver_ATAC-seq/envs/hmmratac.yml"
    shell:
        """
        source activate macs3
        macs3 hmmratac --bam {input.bam} --outdir {config[peak_dir]} --name {wildcards.sample} --save-digested --save-states
        """
#-i {input.index} -g {input.genome} --window {config[hmmratac][window]}

rule filter_peaks:
    resources:
        nodes=1,
        task=1,
        mem='10G',
        time='04:00:00',
    input:
        os.path.join(config['peak_dir'],'{sample}_peaks.gappedPeak')
    output:
        os.path.join(config['peak_dir'],'{sample}.filteredPeaks.gappedPeak')
    shell:
        """
        awk -v OFS="\t" '$13>=10 {{print}}' {input} > {output}
        """

rule filter_summits:
    resources:
        nodes=1,
        task=1,
        mem='10G',
        time='04:00:00',
    input:
        os.path.join(config['peak_dir'],'{sample}_summits.bed')
    output:
        os.path.join(config['peak_dir'],'{sample}.filteredSummits.bed')
    shell:
        """
        awk -v OFS="\t" '$5>=10 {{print}}' {input} > {output}
        """









































