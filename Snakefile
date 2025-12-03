# Snakefile
import pandas as pd
configfile: "config.yaml"


meta = pd.read_csv(config["metadata"], sep="\t")
SAMPLES = meta["sample"].tolist()
META_DICT = meta.set_index("sample").to_dict(orient="index")

OUTDIR = config["outdir"]

rule all:
    input:
        os.path.join(OUTDIR, "featurecounts", "featureCounts.txt"),
        os.path.join(OUTDIR, "plots", "counts_summary.pdf")


rule fastp:
    input:
        R1=lambda wildcards: META_DICT[wildcards.sample]["R1"],
        R2=lambda wildcards: META_DICT[wildcards.sample]["R2"]
    output:
        r1 = temp(os.path.join(OUTDIR, "fastp", "{sample}_R1.clean.fastq.gz")),
        r2 = temp(os.path.join(OUTDIR, "fastp", "{sample}_R2.clean.fastq.gz")),
        json = os.path.join(OUTDIR, "fastp", "{sample}.json"),
        html = os.path.join(OUTDIR, "fastp", "{sample}.html")
    threads: config["fastp"]["threads"]
    params:
        cut_front = "--cut_front" if config["fastp"]["cut_front"] else "",
        cut_front_window_size = f"--cut_front_window_size {config['fastp']['cut_front_window_size']}",
        cut_front_mean_quality = f"--cut_front_mean_quality {config['fastp']['cut_front_mean_quality']}",
        cut_tail = "--cut_tail" if config["fastp"]["cut_tail"] else "",
        cut_tail_window_size = f"--cut_tail_window_size {config['fastp']['cut_tail_window_size']}",
        cut_tail_mean_quality = f"--cut_tail_mean_quality {config['fastp']['cut_tail_mean_quality']}",
        overlap_len_require = f"--overlap_len_require {config['fastp']['overlap_len_require']}",
        length_required = f"--length_required {config['fastp']['length_required']}"
    shell:
        r"""
        fastp \
            -i {input.R1} -I {input.R2} \
            -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} \
            -w {threads} \
            {params.cut_front} {params.cut_front_window_size} {params.cut_front_mean_quality} \
            {params.cut_tail} {params.cut_tail_window_size} {params.cut_tail_mean_quality} \
            {params.overlap_len_require} {params.length_required}
        """

rule bowtie2_align:
    input:
        r1 = os.path.join(OUTDIR, "fastp", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(OUTDIR, "fastp", "{sample}_R2.clean.fastq.gz")
    output:
        bam = temp(os.path.join(OUTDIR, "align", "{sample}.unsorted.bam"))
    threads: config["threads"]["align"]
    params:
        index = config["bowtie2_index"]
    shell:
        r"""
        bowtie2 --threads {threads} --end-to-end --very-sensitive \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
        | samtools view -@ {threads} -b -F 4 -o {output.bam} -
        """

############################
#      SAMTOOLS SORT      #
############################
rule samtools_sort:
    input:
        os.path.join(OUTDIR, "align", "{sample}.unsorted.bam")
    output:
        os.path.join(OUTDIR, "align", "{sample}.sorted.bam")
    threads: config["threads"]["samtools"]
    shell:
        r"""
        samtools sort -@ {threads} -o {output} {input}
        """

rule samtools_index:
    input:
        os.path.join(OUTDIR, "align", "{sample}.sorted.bam")
    output:
        os.path.join(OUTDIR, "align", "{sample}.sorted.bam.bai")
    shell:
        r"""
        samtools index {input}
        """
rule featurecounts:
    input:
        bam = expand(os.path.join(OUTDIR, "align", "{sample}.sorted.bam"), sample=SAMPLES),
        bai = expand(os.path.join(OUTDIR, "align", "{sample}.sorted.bam.bai"), sample=SAMPLES)
    output:
        os.path.join(OUTDIR, "featurecounts", "featureCounts.txt")
    threads: config["threads"]["featurecounts"]
    params:
        gtf = config["gtf"],
        outdir = os.path.join(OUTDIR, "featurecounts")
    shell:
        r"""
        mkdir -p {params.outdir}
        featureCounts \
            -T {threads} \
            -a {params.gtf} \
            -o {output} \
            -p -B -t exon -g gene_id \
            {input.bam}
        """

rule plot_counts:
    input:
        os.path.join(OUTDIR, "featurecounts", "featureCounts.txt")
    output:
        os.path.join(OUTDIR, "plots", "counts_summary.pdf")
    params:
        outdir = os.path.join(OUTDIR, "plots")
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript scripts/plot_counts.R {input} {output}
        """
