# wf-cnv

This repository contains a [Nextflow](https://www.nextflow.io/) workflow for carrying out copy number analysis, using a read depth method implemented by the R package QDNAseq. The input to the workflow is sequence data in FASTQ format, and the output per sample is an HTML report containing chromosome copy summary, ideoplot, plot of read counts per bin, links to genes in detected CNVs, and QC data. The workflow also produces read statistics, a BAM alignment file, BED files of both raw and normalised read counts, and a VCF file.

Please note, currently CNV calling is restricted to human genome builds hg19 and hg38. For more information about the workflow, please see [this](https://labs.epi2me.io/copy-number-calling/) EPI2ME labs blog post.
## Introduction

The workflow takes FASTQ sequence data, aligns to a reference genome, and uses the R package QDNAseq to call copy number aberrations.

Best practices for human copy number calling are actively being investigated by the ONT applications team, and this workflow puts some of that work into something that can be easily used by our community.

wf-cnv also utilises our new reporting and plotting package [ezcharts](https://github.com/epi2me-labs/ezcharts). This uses python [dominate](https://github.com/Knio/dominate) and an apache [echart](https://echarts.apache.org/en/index.html) API to allow us to make modern, responsive layouts and plots with relative ease.
## Quickstart

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and
software resources, and as such Nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or [Singularity](https://sylabs.io/singularity/) to provide isolation of the required software. Both methods are automated out-of-the-box provided either Docker or Singularity is installed.


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

The FASTQs for three test samples are available [here](https://github.com/epi2me-labs/wf-cnv/tree/master/test_data/fastq) and can be used with the the accompanying sample sheet from [here](https://github.com/epi2me-labs/wf-cnv/blob/master/test_data/sample_sheet.csv).

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
* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/products/docker-desktop)
* [Singularity](https://sylabs.io/singularity/)

**Reference**

Scheinin I, Sie D, Bengtsson H, van de Wiel MA, Olshen AB, van Thuijl HF, van Essen HF, Eijk PP, Rustenburg F, Meijer GA, Reijneveld JC, Wesseling P, Pinkel D, Albertson DG, Ylstra B. DNA copy number analysis of fresh and formalin-fixed specimens by shallow whole-genome sequencing with identification and exclusion of problematic regions in the genome assembly. Genome Res. 2014 Dec;24(12):2022-32. doi: 10.1101/gr.175141.114. Epub 2014 Sep 18. [PMCID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4248318/)
