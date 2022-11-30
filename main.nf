#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'

process concatenateReads {
    // concatenate fastq and fastq.gz in a dir
    label params.process_label
    cpus 1
    input:
        tuple path(directory), val(meta)
        val(genome)
    output:
      tuple val(meta.sample_id), val(meta.type), path("${meta.sample_id}.fastq.gz")
      tuple val(meta.sample_id), path("${meta.sample_id}.stats"), emit: fastqstats
    shell:
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} > ${meta.sample_id}.fastq
    gzip ${meta.sample_id}.fastq
    """
}

process alignment {
  label params.process_label
  cpus params.map_threads
  publishDir "${params.out_dir}/BAM", mode: 'copy', pattern: "*"
  input:
    tuple val(sample_id), val(type), file(fastq)
    file reference

  output:
    tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

  """
  minimap2 -t ${task.cpus} -ax map-ont ${reference} ${fastq} | samtools sort -o ${sample_id}.bam
  samtools index ${sample_id}.bam
  """
}


process callCNV {
  label params.process_label
  cpus 1
  publishDir "${params.out_dir}/qdna_seq", mode: 'copy', pattern: "*"
  input:
    tuple val(sample_id), val(type), file(bam), file(bai)

  output:
    tuple val(sample_id), val(type), path("${sample_id}_combined.bed"), path("${sample_id}*"), path("${sample_id}_noise_plot.png"), path("${sample_id}_isobar_plot.png")

  script:
  """
  run_qdnaseq.r --bam ${bam} --out_prefix ${sample_id} --binsize ${params.bin_size} --reference ${params.genome}
  cut -f5 ${sample_id}_calls.bed | paste ${sample_id}_bins.bed - > ${sample_id}_combined.bed
  """
}

process checkFASTA {
    label params.process_label
    cpus 1
    input:
      file(reference)
    output:
      path("${reference}_genome.txt")
      path("output.txt"), emit: genome_label optional true

    script:
    """
    faSize -detailed -tab ${reference} > ${reference}_genome.txt
    check_fasta.py --chr_counts ${reference}_genome.txt --genome ${params.genome} -o output.txt
    """
}

//

process makeReport {
  label params.process_label
  cpus 1
  input:
    tuple val(sample_id), path(read_stats), val(type), path(cnv_calls), val(cnv_files), path(noise_plot), path(isobar_plot)
    path "versions/*"
    path "params.json"
  output:
    path("${sample_id}_wf-cnv-report.html")
  script:
  """
  cnv_plot.py \
    -q ${cnv_calls} \
    -o ${sample_id}_wf-cnv-report.html \
    --read_stats ${read_stats}\
    --params params.json \
    --versions versions \
    --bin_size ${params.bin_size} \
    --genome ${params.genome} \
    --sample_id ${sample_id} \
    --noise_plot ${noise_plot} \
    --isobar_plot ${isobar_plot}
  """

}

process getVersions {
    label params.process_label
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    R --version | grep -w R | grep version | cut -f3 -d" " | sed 's/^/R,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    R --slave -e 'packageVersion("QDNAseq")' | cut -d\\' -f2 | sed 's/^/QDNAseq,/' >> versions.txt
    """
}


process getParams {
    label params.process_label
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label params.process_label
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}


// workflow module
workflow pipeline {
    take:
        reads
        reference
    main:

        genome = checkFASTA(reference)

        genome_match_channel = genome.genome_label.ifEmpty{exit 1, log.error('Reference FASTA and selected genome do not match')}

        software_versions = getVersions()
        workflow_params = getParams()

        sample_fastqs = concatenateReads(reads, genome.genome_label)

        alignment = alignment(sample_fastqs[0], reference)

        cnvs = callCNV(alignment)

        join_ch = sample_fastqs.fastqstats.join(cnvs.groupTuple())

        report = makeReport(join_ch, software_versions.collect(), workflow_params)

    emit:
        results = report.concat(workflow_params)
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

valid_bin_sizes = [1, 5, 10, 15, 30, 50, 100, 500, 1000]
valid_genomes = ["hg19", "hg38"]


workflow {

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

     if (workflow.profile == "conda") {
       log.error "conda is not supported by this workflow"
       exit 1
     }

    if (!valid_bin_sizes.any { it == params.bin_size}) {
      log.error "`--bin-size` should be one of: $valid_bin_sizes"
      exit 1
    }

    if (!valid_genomes.any { it == params.genome}) {
     log.error "`--genome` should be one of: $valid_genomes"
     exit 1
   }

    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet])

    if (params.reference) {
      reference_fasta = file(params.reference, checkIfExists: true)
    } else {
      exit 1, 'Fasta file not specified!'
    }

    pipeline(samples, reference_fasta)
    output(pipeline.out.results)
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
