# Workflow template

This repository contains a [nextflow](https://www.nextflow.io/) workflow for carrying out copy number analysis.

> This workflow is not intended to be used by end users.
## Introduction

The workflow takes FASTQ sequence data, aligns to a reference genome, and uses the R package QDNAseq to call copy number aberrations.

Best practices for human copy number calling are actively being investigated by the ONT applications team, and this workflow puts some of that work into something that can be easily used by our community.

wf-cnv also utilises our new reporting and plotting package [ezcharts](https://github.com/epi2me-labs/ezcharts). This uses python [dominate](https://github.com/Knio/dominate) and an apache [echart](https://echarts.apache.org/en/index.html) api to allow us to make modern, responsive layouts and plots with relative ease.

With the release of wf-cnv we also include a new ideogram plotting component for ezharts.
## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-cnv --help
```

to see the options for the workflow.

Example command:

```
nextflow run epi2me-labs/wf-cnv --fastq <PATH_TO_FASTQS> --fasta <PATH_TO_REFERENCE> --genome <hg19|hg38> --bin_size <BIN_SIZE>
```

The FASTQs for three test samples are available [here](https://github.com/epi2me-labs/wf-cnv/test_data/fastq) and can be used with the the accompanying sample sheet from [here](https://github.com/epi2me-labs/wf-cnv/test_data/sample_sheet.csv).

Example command with test data:

```
nextflow run epi2me-labs/wf-cnv --fastq <PATH_TO_DOWNLOADED_FASTQ> --sample_sheet <PATH_TO_DOWNLOADED_SAMPLE_SHEET> --fasta /path/to/hg38.fa.gz --genome hg38 --bin_size 500
```

**Workflow outputs**

The primary outputs of the workflow include, per sample:

* `<SAMPLE_NAME>_wf-cnv-report.html`: HTML CNV report containing chromosome copy summary, ideoplot, plot of read counts per bin, links to genes in detected CNVs, and QC data: read length histogram, noise plot (noise as a function of sequence depth) and isobar plot (median read counts per bin shown as a function of GC content and mappability)
* `<SAMPLE_NAME>.stats`: Read stats
* `BAM/<SAMPLE_NAME>.bam`: Alignment of reads to reference
* `BAM/<SAMPLE_NAME>.bam.bai`: BAM index
* `qdna_seq/<SAMPLE_NAME>_plots.pdf`: QDNAseq-generated plots
* `qdna_seq/<SAMPLE_NAME>_raw_bins.bed`: BED file of raw read counts per bin
* `qdna_seq/<SAMPLE_NAME>_bins.bed`: BED file of corrected, normalised, and smoothed read counts per bin
* `qdna_seq/<SAMPLE_NAME>_calls.vcf`: VCF file of CNV calls
## Useful links


* [QDNAseq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html)
* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
