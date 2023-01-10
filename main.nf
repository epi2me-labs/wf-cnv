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
include { sample_ingress } from './lib/ingress'

process concatenateReads {
    // concatenate fastq and fastq.gz in a dir
    label params.process_label
    cpus 1
    input:
        tuple path(directory), val(meta)
        val(genome)
    output:
      tuple val(meta.sample_id), val(meta.type), path("${meta.sample_id}.fastq.gz")
      tuple val(meta.sample_id), path("${meta.sample_id}.stats"), emit: stats
    shell:
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} > ${meta.sample_id}.fastq
    gzip ${meta.sample_id}.fastq
    """
}

process bamStats {
    label params.process_label
    cpus 1
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
    output:
        tuple val(sample_id), path("${sample_id}.stats"), emit: stats
    shell:
    """
    bamstats --threads 3 ${sample_id}.bam > ${sample_id}.stats
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
        tuple val(sample_id), val(type), path("${sample_id}.bam")
    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${reference} ${fastq} | samtools sort -o ${sample_id}.bam
    """
}

process samtoolsIndex {
    label params.process_label
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam")
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
    script:
    """
    samtools index ${sample_id}.bam
    """
}

process callCNV {
    label params.process_label
    cpus 1
    publishDir "${params.out_dir}/qdna_seq", mode: 'copy', pattern: "*"
    input:
        tuple val(sample_id), val(type), file(bam), file(bai)
        val(genome_build)
    output:
    tuple val(sample_id), val(type), path("${sample_id}_combined.bed"), path("${sample_id}*"), path("${sample_id}_noise_plot.png"), path("${sample_id}_isobar_plot.png")

    script:
    """
    run_qdnaseq.r --bam ${bam} --out_prefix ${sample_id} --binsize ${params.bin_size} --reference ${genome_build}
    cut -f5 ${sample_id}_calls.bed | paste ${sample_id}_bins.bed - > ${sample_id}_combined.bed
    """
}

process getGenome {
    label params.process_label
    cpus 1
    input:
        file(reference)
    output:
        path("${reference}_genome.txt")
        path("output.txt"), emit: genome_label optional true
        env genome_string, emit: genome_build optional true
     script:
    if (params.fastq) 
        """
        faSize -detailed -tab ${reference} > ${reference}_genome.txt
        get_genome.py --chr_counts ${reference}_genome.txt -o output.txt
        genome_string=`cat output.txt`
        """
    else if(params.bam)
        """
        samtools idxstats ${reference} > ${reference}_genome.txt
        get_genome.py --chr_counts ${reference}_genome.txt -o output.txt
        genome_string=`cat output.txt`
        """
}


process makeReport {
    label params.process_label
    cpus 1
    input:
        tuple val(sample_id), path(read_stats), val(type), path(cnv_calls), val(cnv_files), path(noise_plot), path(isobar_plot)
        path "versions/*"
        path "params.json"
        val genome_build
    output:
        path("wf-cnv-report.html")

    script:
    """
    cnv_plot.py \
        -q ${cnv_calls} \
        -o wf-cnv-report.html \
        --read_stats ${read_stats}\
        --params params.json \
        --versions versions \
        --bin_size ${params.bin_size} \
        --genome ${genome_build} \
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

        software_versions = getVersions()
        workflow_params = getParams()

        if (params.fastq) {

            genome = getGenome(reference)
            //empty genome.genome_build means chromosome sizes didn't match to either hg19 or hg38, exit
            genome_match_channel = genome.genome_build.ifEmpty{exit 1, log.error('The genome build detected is not compatible with this workflow. Plesae check the reference sequence and try again.')}
            sample_fastqs = concatenateReads(reads, genome.genome_label) 
            stats = sample_fastqs.stats
            alignment = alignment(sample_fastqs[0], reference) 
            bam_for_cnv = samtoolsIndex(alignment)

        } else if (params.bam) {

            bam_file = Channel.fromPath(params.bam) 
            bam_idx_fp = params.bam + ".bai" 
            sample_details = reads.map{it -> tuple(it[1].sample_id, it[1].type)}.flatten()
            
            if(file(bam_idx_fp).exists()) {

                bam_idx_file = Channel.fromPath(bam_idx_fp)
                bam_for_cnv = sample_details.concat(bam_file, bam_idx_file).collate(4)
                
            } else {
                
                bam_channel = sample_details.concat(bam_file).collate(3)
                bam_for_cnv = samtoolsIndex(bam_channel) 

            }

            genome = getGenome(bam_file)
            genome_match_channel = genome.genome_build.ifEmpty{exit 1, log.error('The genome build detected is not compatible with this workflow. Plesae check the reference sequence and try again.')}
            stats = bamStats(bam_for_cnv)
        }
        
        cnvs = callCNV(bam_for_cnv, genome.genome_build) 
        join_ch = stats.join(cnvs.groupTuple())
        report = makeReport(join_ch, software_versions.collect(), workflow_params, genome.genome_build)


    emit:
        results = report.concat(workflow_params)
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

valid_bin_sizes = [1, 5, 10, 15, 30, 50, 100, 500, 1000]

workflow {

    if (params.bam && params.fastq) {
        log.error "Please choose only one of --fastq or --bam"
        exit 1
    }

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

    if (params.fastq != null) {
        
        // Reference file required when using FASTQ
        // reference is determined from BAM header when using BAM
        if (params.reference) {

            reference_fasta = file(params.reference, checkIfExists: true)

        } else {
            
            log.error "Reference file not specified!"
            exit 1

        }

        samples = fastq_ingress([
          "input": params.fastq,
          "sample": params.sample,
          "sample_sheet": params.sample_sheet])

    
    } else if (params.bam != null) {

        reference_fasta = null

        samples = sample_ingress([
          "input": params.bam,
          "sample": params.sample,
          "sample_sheet": params.sample_sheet,
          "input_type": "BAM"])

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
