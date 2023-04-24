import os
import pandas as pd

#configfile: "config/Liver_male_female_atac_seq.yaml"
configfile: "config/All_tissue_male_female_atac_seq.yaml"
READS = ['1', '2']
sample_metadata = pd.read_csv(config['sample_metadata'])

rule all:
    input:
        #expand(os.path.join(config['fastqc_dir'],"{sample}_{read}_fastqc.html"),
        #    sample=config['SAMPLE'],
        #    read=READS),
        expand(os.path.join(config['mapping_dir'],'{sample}.sam'),
            sample=config['SAMPLE'])
        #expand(os.path.join(config['mapping_dir'],'{sample}_coverage.bw'),
        #    sample=config['SAMPLE'])
        #expand(os.path.join(config['mapping_dir'],'{tissue}_coverage.tsv'),
        #    tissue=config['TISSUE'])
        #expand(os.path.join(config['peak_dir'],'{sample}.filteredPeaks.gappedPeak'),
        #    sample=config['SAMPLE']),
        #expand(os.path.join(config['peak_dir'],'{sample}.filteredSummits.bed'),
        #    sample=config['SAMPLE'])

rule fastqc:
    input:
        os.path.join(config['input_dir'],"{sample}_{read}.fastq")
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
        [os.path.join(config['input_dir'],"{{sample}}_{read}.fastq").format(read=read) for read in READS]
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
    threads: 24
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
    threads: 4
    shell:
        "bamCoverage -b {input} -o {output} --normalizeUsing CPM -p {threads}"

rule combine_bigwigs:
    input:
        lambda wildcards: [os.path.join(config['mapping_dir'],f'{s}_coverage.bw') for s in sample_metadata.loc[sample_metadata.TISSUE==wildcards.tissue,'Run']]
    params:
        labels = lambda wildcards: list(sample_metadata.loc[sample_metadata.TISSUE==wildcards.tissue,'Sample Name'])
    output:
        npz=os.path.join(config['mapping_dir'],"{tissue}_coverage.npz"),
        tsv=os.path.join(config['mapping_dir'],"{tissue}_coverage.tsv")
    threads: 4
    shell:
        "multiBigwigSummary bins -b {input} -p {threads} -bs 50 -o {output.npz} --outRawCounts {output.tsv} --labels {params.labels}"

rule peak_calling:
    input:
        bam = os.path.join(config['mapping_dir'],'{sample}_sorted.bam'),
        index = os.path.join(config['mapping_dir'],'{sample}_sorted.bam.bai'),
        genome = os.path.join(config['mapping_dir'],'{sample}_genome.info')
    output:
        os.path.join(config['peak_dir'],'{sample}_peaks.gappedPeak'),
        os.path.join(config['peak_dir'],'{sample}_summits.bed')
    shell:
        "HMMRATAC -b {input.bam} -i {input.index} -g {input.genome} --window {config[hmmratac][window]}"

rule filter_peaks:
    input:
        os.path.join(config['peak_dir'],'{sample}_peaks.gappedPeak')
    output:
        os.path.join(config['peak_dir'],'{sample}.filteredPeaks.gappedPeak')
    shell:
        """
        awk -v OFS="\t" '$13>=10 {{print}}' {input} > {output}
        """

rule filter_summits:
    input:
        os.path.join(config['peak_dir'],'{sample}_summits.bed')
    output:
        os.path.join(config['peak_dir'],'{sample}.filteredSummits.bed')
    shell:
        """
        awk -v OFS="\t" '$5>=10 {{print}}' {input} > {output}
        """



