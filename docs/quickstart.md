## Quickstart

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and
software resources, and as such Nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or [Singularity](https://sylabs.io/singularity/) to provide isolation of the required software. Both methods are automated out-of-the-box provided either Docker or Singularity is installed.


It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-cnv --help
```

to see the options for the workflow.

Example command (BAM):

```
nextflow run epi2me-labs/wf-cnv --bam <PATH_TO_BAM> --bin_size <BIN_SIZE>
```

Example command (FASTQ):

```
nextflow run epi2me-labs/wf-cnv --fastq <PATH_TO_FASTQS> --reference <PATH_TO_REFERENCE> --bin_size <BIN_SIZE>
```

The FASTQs for three test samples are available [here](https://github.com/epi2me-labs/wf-cnv/tree/master/test_data/fastq) and can be used with the the accompanying sample sheet from [here](https://github.com/epi2me-labs/wf-cnv/blob/master/test_data/sample_sheet.csv).

Example command with test data:

```
nextflow run epi2me-labs/wf-cnv --fastq <PATH_TO_DOWNLOADED_FASTQ> --sample_sheet <PATH_TO_DOWNLOADED_SAMPLE_SHEET> --reference /path/to/hg38.fa.gz --bin_size 500
```

**Workflow outputs**

The primary outputs of the workflow include, per sample:

* `<SAMPLE_NAME>_wf-cnv-report.html`: HTML CNV report containing chromosome copy summary, ideoplot, plot of read counts per bin, links to genes in detected CNVs, and QC data: read length histogram, noise plot (noise as a function of sequence depth) and isobar plot (median read counts per bin shown as a function of GC content and mappability)
* `<SAMPLE_NAME>.stats`: Read stats
* `BAM/<SAMPLE_NAME>.bam`: Alignment of reads to reference (FASTQ input)
* `BAM/<SAMPLE_NAME>.bam.bai`: BAM index (FASTQ input)
* `qdna_seq/<SAMPLE_NAME>_plots.pdf`: QDNAseq-generated plots
* `qdna_seq/<SAMPLE_NAME>_raw_bins.bed`: BED file of raw read counts per bin
* `qdna_seq/<SAMPLE_NAME>_bins.bed`: BED file of corrected, normalised, and smoothed read counts per bin
* `qdna_seq/<SAMPLE_NAME>_calls.vcf`: VCF file of CNV calls
