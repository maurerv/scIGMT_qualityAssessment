# Assess quality of bisulfite and rna sequencing data

## Overview

This pipeline uses the following tools to assess data quality:

1. [FastQC](https://github.com/s-andrews/FastQC), to assess fastq data quality
2. [samtools](https://github.com/samtools/) for alignment statistics
3. [mosdepth](https://github.com/brentp/mosdepth) to calculate alignment depth
4. [picard](https://github.com/broadinstitute/picard) to determine insert size distributions
5. [methyldackel](https://github.com/dpryan79/MethylDackel) to compute bisulfite conversion rates
6. [MultiQC](https://multiqc.info/) to compile everything into a neat report


## Configuration

All necessary tools can be installed using conda and the provided environment file as per:
```
conda env create -f env.yml
```
The environment specification file does not control the versions of the utilized tools, so its up to the user to keep versions consistent across datasets.

All pipeline steps can also be run within singularity containers, for which [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) needs to be installed and configured.


## Running the pipeline

Adapting the pipeline to your use cases requires at least three changes to the config.yaml.


1. ```input``` is the absolute path to a semicolon separated sample sheet without headers and the following format:
    | Modality  | Path | ID |
    | ----- | ---- | ---- |
    | Methylome or Transcriptome  | Absolute path to file | Unique sample ID |

2. ```OUTDIR``` is the absolute path to directory to which the results will be written,
3. ```REFERENCE``` is the absolute path to a reference genome fasta file.

```FASTQ_EXTENSIONS``` and ```BAM_EXTENSIONS``` are optional files that specify the range of file suffixes that mark a given file as FASTQ or BAM.


## Result

After running the pipeline, the output directory will have the following structure
```
OUTDIR
├── Methylome/            # Contains per sample quality statistics and report for Methylome
├── Transcriptome/        # Contains per sample quality statistics and report for Transcriptome
```
