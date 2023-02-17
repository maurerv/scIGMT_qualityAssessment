from os import extsep
from os.path import join, realpath, exists, splitext, basename

configfile : "config.yaml"

OUTDIR           = config["OUTDIR"]
FASTQ_EXTENSIONS = tuple(config["FASTQ_EXTENSIONS"])
BAM_EXTENSIONS   = tuple(config["BAM_EXTENSIONS"])

# Read in sample sheet, check if files exist and determine whether they are bam/fastq
SAMPLE_SHEET = {"Modality" : [], "Path" : [], "Type" : [], "ID" : []}
with open(config["input"], "r", encoding = "utf-8") as ifile:
    for index, line in enumerate(ifile):
        line = line.split(";")
        if len(line) != 3:
            raise ValueError(f"Line {index} does not follow the expected format.")
        path, modality, sample_id = [element.strip() for element in line]
        if not exists(path):
            raise FileNotFoundError(f"Line {index} contains non existing file.")
        if path.endswith(FASTQ_EXTENSIONS):
            SAMPLE_SHEET["Type"].append("fastq")
        elif path.endswith(BAM_EXTENSIONS):
            SAMPLE_SHEET["Type"].append("bam")
        else:
            raise NotImplementedError(f"Filetype in line {index} is not supported.")
        if modality not in ("Methylome", "Transcriptome"):
            raise ValueError("Only 'Methylome''Transcriptome' are accepted modalities")
        SAMPLE_SHEET["Modality"].append(modality)
        SAMPLE_SHEET["ID"].append(sample_id)
        SAMPLE_SHEET["Path"].append(path)

indices = [i for i, v in enumerate(SAMPLE_SHEET["Type"]) if v == "fastq"]
fastqs  = [SAMPLE_SHEET["Path"][i] for i in indices]
fastq_modality  = [SAMPLE_SHEET["Modality"][i] for i in indices]
fastq_samples = [basename(fastq).split(extsep)[0] for fastq in fastqs]
fastq_in = dict(zip(zip(fastq_samples, fastq_modality), fastqs))

indices = [i for i, v in enumerate(SAMPLE_SHEET["Type"]) if v == "bam"]
bams  = [SAMPLE_SHEET["Path"][i] for i in indices]
bam_modality  = [SAMPLE_SHEET["Modality"][i] for i in indices]
bam_samples = [splitext(basename(bam))[0] for bam in bams]
bam_in = dict(zip(zip(bam_samples, bam_modality), bams))

bams_bs_idx      = [i for i in indices if SAMPLE_SHEET["Modality"][i] == "Methylome"]
bams_bs          = [SAMPLE_SHEET["Path"][i] for i in bams_bs_idx]
bam_modality_bs  = [SAMPLE_SHEET["Modality"][i] for i in bams_bs_idx]
bam_samples_bs   = [splitext(basename(bam))[0] for bam in bams_bs]

modalities = set(SAMPLE_SHEET["Modality"])

def find_fastqs(wildcards):
    return fastq_in[(wildcards.sample, wildcards.modality)]

def find_bams(wildcards):
    return bam_in[(wildcards.sample, wildcards.modality)]


rule all:
    input:
        expand("{outdir}/{modality}/multiqc_report.html",
            zip, outdir = [OUTDIR] * len(modalities), modality = modalities),
        expand("{outdir}/{modality}/{sample}.conversionRate.txt",
            zip, outdir = [OUTDIR] * len(bams_bs), modality = bam_modality_bs,
            sample = bam_samples_bs
        )

rule fastqc:
    input:
        find_fastqs,
    output:
        join(OUTDIR, "{modality}", "{sample}_fastqc.zip"),
    singularity:
        "docker://pegi3s/fastqc"
    resources:
        avg_mem  = lambda wildcards, attempt: 800 * attempt,
        mem_mb   = lambda wildcards, attempt: 1000 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    threads: 1
    params:
        outdir = OUTDIR
    shell:
        "fastqc {input} -o {params.outdir}/{wildcards.modality}"

rule samtools_stats:
    input:
        find_bams,
    output:
        join(OUTDIR, "{modality}", "{sample}.samtoolsStats.txt")
    singularity:
        "docker://staphb/samtools"
    resources:
        avg_mem  = lambda wildcards, attempt: 400 * attempt,
        mem_mb   = lambda wildcards, attempt: 600 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    threads: 1
    shell:
        "samtools stats {input} > {output}"

rule bam_coverage:
    input:
        find_bams,
    output:
        join(OUTDIR, "{modality}", "{sample}.mosdepth.global.dist.txt")
   singularity:
        "docker://quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
    resources:
        avg_mem  = lambda wildcards, attempt: 2000 * attempt,
        mem_mb   = lambda wildcards, attempt: 2400 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    threads: 1
    params:
        outdir = OUTDIR
    shell:
        "mosdepth -n {params.outdir}/{wildcards.modality}/{wildcards.sample} {input}"

rule bisulfite_qc:
    input:
        find_bams,
    output:
        convRate = join(OUTDIR, "{modality}", "{sample}.conversionRate.txt"),
        cpgStats = join(OUTDIR, "{modality}", "{sample}.CpGCoverage.txt")
    singularity:
        "docker://nfcore/methylseq"
    params:
        genome_reference = config["REFERENCE"]
    resources:
        avg_mem  = lambda wildcards, attempt: 1000 * attempt,
        mem_mb   = lambda wildcards, attempt: 1200 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    threads: 4
    shell:"""
        MethylDackel extract \
            --CHH \
            -o {output.convRate} \
            -@ {threads} \
            {params.genome_reference} \
            {input}

        awk '
          BEGIN{{FS=OFS="\t"}}
          NR>1{{
            rate[$1]+=$4
            meth[$1]+=$5
            umeth[$1]+=$6
            cov[$1]+=1
          }}
          END{{
            for (i in meth)
              print i,cov[i],meth[i],umeth[i],1-(meth[i]/(meth[i]+umeth[i])),rate[i]/cov[i]
          }}' {output.convRate}_CpG.bedGraph > {output.cpgStats}
        rm {output.convRate}_CpG.bedGraph

        awk '
          BEGIN{{FS=OFS="\t"}}
          NR>1{{
            rate[$1]+=$4
            meth[$1]+=$5
            umeth[$1]+=$6
            cov[$1]+=1
          }}
          END{{
            for (i in meth)
              print i,cov[i],meth[i],umeth[i],1-(meth[i]/(meth[i]+umeth[i])),rate[i]/cov[i]
          }}' {output.convRate}_CHH.bedGraph > {output.convRate}
        rm {output.convRate}_CHH.bedGraph

        """

rule insert_size:
    input:
        find_bams,
    output:
        join(OUTDIR, "{modality}", "{sample}.insertSize.txt")
    singularity:
        "docker://dquz/picard:latest"
    resources:
        avg_mem  = lambda wildcards, attempt: 800 * attempt,
        mem_mb   = lambda wildcards, attempt: 1000 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    threads: 1
    shell:"""
        picard CollectInsertSizeMetrics \
            -I {input} \
            -O {output} \
            -H /dev/null \
            -M 0.5
        """

rule multiqc:
    input:
        expand("{outdir}/{modality}/{sample}_fastqc.zip",
            zip, outdir = [OUTDIR] * len(fastqs), modality = fastq_modality,
            sample = fastq_samples
        ),
        expand("{outdir}/{modality}/{sample}.samtoolsStats.txt",
            zip, outdir = [OUTDIR] * len(bams), modality = bam_modality,
            sample = bam_samples
        ),
        expand("{outdir}/{modality}/{sample}.mosdepth.global.dist.txt",
            zip, outdir = [OUTDIR] * len(bams), modality = bam_modality,
            sample = bam_samples
        ),
        expand("{outdir}/{modality}/{sample}.insertSize.txt",
            zip, outdir = [OUTDIR] * len(bams), modality = bam_modality,
            sample = bam_samples
        )
    output:
        join(OUTDIR, "{modality}", "multiqc_report.html")
    singularity:
        "docker://ewels/multiqc"
    resources:
        avg_mem  = lambda wildcards, attempt: 5000 * attempt,
        mem_mb   = lambda wildcards, attempt: 6000 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    threads: 1
    params:
        idir = join(OUTDIR, "{modality}")
    shell:"""
        multiqc -f -o {params.idir} {params.idir}
        """
