# wf-cnv

This repository contains a [Nextflow](https://www.nextflow.io/) workflow for carrying out copy number analysis, using a read depth method implemented by the R package QDNAseq. The input to the workflow is sequence data in FASTQ format, and the output per sample is an HTML report containing chromosome copy summary, ideoplot, plot of read counts per bin, links to genes in detected CNVs, and QC data. The workflow also produces read statistics, a BAM alignment file, BED files of both raw and normalised read counts, and a VCF file.

Please note, currently CNV calling is restricted to human genome builds hg19 and hg38. For more information about the workflow, please see [this](https://labs.epi2me.io/copy-number-calling/) EPI2ME labs blog post.
