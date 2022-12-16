## Introduction

The workflow takes BAM or FASTQ data, aligns to a reference genome (if FASTQ files are supplied), and uses the R package QDNAseq to call copy number aberrations.

Best practices for human copy number calling are actively being investigated by the ONT applications team, and this workflow puts some of that work into something that can be easily used by our community.

wf-cnv also utilises our new reporting and plotting package [ezcharts](https://github.com/epi2me-labs/ezcharts). This uses Python [dominate](https://github.com/Knio/dominate) and an Apache [echart](https://echarts.apache.org/en/index.html) API to allow us to make modern, responsive layouts and plots with relative ease.
