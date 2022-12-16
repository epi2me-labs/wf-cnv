import ArgumentParser

process handleSingleFile {
    label params.process_label
    cpus 1
    input:
        file reads
    output:
        path "$reads.simpleName"
    script:
        def name = reads.simpleName
        def reads_dir = 'reads_dir'
    """
    mkdir $name
    mv $reads $name
    """
}


process checkSampleSheet {
    label params.process_label
    cpus 1
    input:
        file "sample_sheet.txt"
    output:
        file "samples.txt"
    """
    check_sample_sheet.py sample_sheet.txt samples.txt
    """
}

/**
 * Compare number of samples in samplesheet
 * with the number of barcoded dirs found and
 * print warnings
 *
 *
 * @param number of samples in sample sheet
 * @param number of barcoded directories
 * @return null
 */

def compareSampleSheet(int sample_sheet_count, int valid_dir_count)
{

    if (sample_sheet_count != valid_dir_count) {
        log.warn """The number of samplesheet entries ({}) does not match the number of barcoded directories ({})""", sample_sheet_count, valid_dir_count
    }

}

/**
 * Take an input file and sample name to return a channel with
 * a single named sample.
 *
 *
 * @param input_file Single file
 * @param sample_name Name to give the sample
 * @return Channel of tuples (path, map(sample_id, type, barcode))
 */
def handle_single_file(input_file, sample_name)
{
    singleFile = Channel.fromPath(input_file)
    sample = handleSingleFile(singleFile)
    return sample.map { it -> tuple(it, create_metamap([sample_id:sample_name ?: it.simpleName])) }

}


/**
 * Find data using various globs. Wrapper around Nextflow `file`
 * method.
 *
 * @param pattern file object corresponding to top level input folder.
 * @param maxdepth maximum depth to traverse
 * @return list of files.
 */

def find_files(pattern, maxdepth, input_type)
{
    files = []
    def extension_map = ["fastq": ["fastq", "fastq.gz", "fq", "fq.gz"], "bam": ["bam"], "ubam": ["ubam"]]
    extensions = extension_map[input_type]
    for (ext in extensions) {
        files += file(pattern.resolve("*.${ext}"), type: 'file', maxdepth: maxdepth)
    }
    return files
}


/**
 * Rework EPI2ME flattened directory structure into standard form
 * files are matched on barcode\d+ and moved into corresponding
 * subdirectories ready for processing.
 *
 * @param input_folder Top-level input directory.
 * @param staging Top-level output_directory.
 * @return A File object representating the staging directory created
 *     under output
 */
def sanitize_sample(input_folder, staging, input_type)
{
    // TODO: this fails if input_folder is an S3 path
    println("Running sanitization.")
    println(" - Moving files: ${input_folder} -> ${staging}")
    staging.mkdirs()
    files = find_files(input_folder.resolve("**"), 1, input_type)
    for (sample in files) {
        fname = sample.getFileName()
        // find barcode
        pattern = ~/barcode\d+/
        matcher = fname =~ pattern
        if (!matcher.find()) {
            // not barcoded - leave alone
            sample.renameTo(staging.resolve(fname))
        } else {
            bc_dir = file(staging.resolve(matcher[0]))
            bc_dir.mkdirs()
            sample.renameTo(staging.resolve("${matcher[0]}/${fname}"))
        }
    }
    println(" - Finished sanitization.")
    return staging
}


/**
 * Take an input directory return the barcode and non barcode
 * sub directories contained within.
 *
 *
 * @param input_directory Top level input folder to locate sub directories
 * @return A list containing sublists of barcode and non_barcode sub directories
 */
def get_subdirectories(input_directory)
{
    barcode_dirs = file(input_directory.resolve("barcode*"), type: 'dir', maxdepth: 1)
    all_dirs = file(input_directory.resolve("*"), type: 'dir', maxdepth: 1)
    non_barcoded = ( all_dirs + barcode_dirs ) - all_dirs.intersect(barcode_dirs)
    return [barcode_dirs, non_barcoded]
}


/**
 * Load a sample sheet into a Nextflow channel to map barcodes
 * to sample names.
 *
 * @param samples CSV file according to MinKNOW sample sheet specification
 * @return A Nextflow Channel of tuples (barcode, sample name, sample type)
 */
def get_sample_sheet(sample_sheet)
{
    println("Checking sample sheet.")
    sample_sheet = file(sample_sheet);
    is_file = sample_sheet.isFile()

    if (!is_file) {
        println('Error: `--samples` is not a file.')
        exit 1
    }

    return checkSampleSheet(sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.barcode,
            row.sample_id,
            row.type ? row.type : 'test_sample')
        }
}


/**
 * Take a list of input directories and return directories which are
 * valid, i.e. contains only .fastq(.gz) files.
 *
 *
 * @param input_dirs List of barcoded directories (barcodeXX,mydir...)
 * @return List of valid directories
 */
def get_valid_directories(input_dirs, input_type)
{
    valid_dirs = []
    no_files_dirs = []
    invalid_files_dirs = []
    for (d in input_dirs) {
        valid = true
        sample_files = find_files(d, 1, input_type)
        all_files = file(d.resolve("*"), type: 'file', maxdepth: 1)
        other_files = ( all_files + sample_files ) - all_files.intersect(sample_files)

        if (other_files) {
            valid = false
            invalid_files_dirs << d
        }
        if (!sample_files) {
            valid = false
            no_files_dirs << d
        }
        if (valid) {
            valid_dirs << d
        }
    }
    def extension_map = ["fastq": ["fastq", "fastq.gz", "fq", "fq.gz"], "bam": ["bam"], "ubam": ["ubam"]]
    extensions = extension_map[input_type]
    if (valid_dirs.size() == 0) {
        error_message = "Error: None of the directories given contain $extensions files."
        println(error_message)
        exit 1
    }
    if (no_files_dirs.size() > 0) {
        println("Warning: Excluding directories not containing $extensions files:")
        for (d in no_files_dirs) {
            println("   - ${d}")
        }
    }
    if (invalid_files_dirs.size() > 0) {
        println("Warning: Excluding directories containing non $extensions files:")
        for (d in invalid_files_dirs) {
            println("   - ${d}")
        }
    }
    return valid_dirs
}


/**
 * Take an input directory and sample name to return a channel
 * with a single named sample.
 *
 *
 * @param input_directory Directory of files
 * @param sample_name Name to give the sample
 * @return Channel of tuples (path, map(sample_id, type, barcode))
 */
def handle_flat_dir(input_directory, sample_name, input_type)
{
    valid_dirs= get_valid_directories([ file(input_directory) ], input_type)
    return Channel.fromPath(valid_dirs)
        .map { it -> tuple(it, create_metamap([sample_id:sample_name ?: it.baseName])) }

}


/**
 * Take a list of barcode directories and a sample sheet to return
 * a channel of named samples.
 *
 *
 * @param barcoded_dirs List of barcoded directories (barcodeXX,...)
 * @param sample_sheet List of tuples mapping barcode to sample name
 *     or a simple string for non-multiplexed data.
 * @param min_barcode Minimum barcode to accept.
 * @param max_barcode Maximum (inclusive) barcode to accept.
 * @return Channel of tuples (path, map(sample_id, type, barcode))
 */
def handle_barcoded_dirs(barcoded_dirs, sample_sheet, min_barcode, max_barcode, input_type)
{
    valid_dirs = get_valid_directories(barcoded_dirs, input_type)
    // link sample names to barcode through sample sheet
    if (!sample_sheet) {
        sample_sheet = Channel
            .fromPath(valid_dirs)
            .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
            .filter { barcode_in_range(it, min_barcode, max_barcode) }
            .map { path -> tuple(path.baseName, path.baseName, 'test_sample')}
    } else {

        // return warning if there is a discrepancy between the samplesheet and barcode dirs

        // unclassfied will never be in the sample_sheet so remove
        non_unclassified = valid_dirs
        non_unclassified -= 'unclassified'

        barcode_dirs_found = non_unclassified.size()

        int count = 0

        // We do this instead of .count() because valid_dirs is a list and
        // sample_sheet is a channel - the channel is only populated after
        // checkSampleSheet is complete and so if you compare without
        // waiting for that then the comparisson fails
        sample_sheet_entries = sample_sheet.subscribe onNext: { count++ }, onComplete: { compareSampleSheet(count,barcode_dirs_found) }

    }

    return Channel
        .fromPath(valid_dirs)
        .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
        .filter { barcode_in_range(it, min_barcode, max_barcode) }
        .map { path -> tuple(path.baseName, path) }
        .join(sample_sheet)
        .map { barcode, path, sample, type -> tuple(path, create_metamap([sample_id:sample, type:type, barcode:barcode])) }
}


/**
 * Determine if a barcode path is within a required numeric range
 *
 * @param path barcoded directory (barcodeXX).
 * @param min_barcode Minimum barcode to accept.
 * @param max_barcode Maximum (inclusive) barcode to accept.
 */
def barcode_in_range(path, min_barcode, max_barcode)
{
    pattern = ~/barcode(\d+)/
    matcher = "${path}" =~ pattern
    value = matcher[0][1].toInteger()
    valid = ((value >= min_barcode) && (value <= max_barcode))
    return valid
}


/**
 * Take a list of non-barcode directories to return a channel
 * of named samples. Samples are named by directory baseName.
 *
 *
 * @param non_barcoded_dirs List of directories (mydir,...)
 * @return Channel of tuples (path, map(sample_id, type, barcode))
 */
def handle_non_barcoded_dirs(non_barcoded_dirs, input_type)
{
    valid_dirs = get_valid_directories(non_barcoded_dirs, input_type)
    return Channel.fromPath(valid_dirs)
        .map { path -> tuple(path, create_metamap([sample_id:path.baseName])) }
}

def create_metamap(Map arguments) {
    def parser = new ArgumentParser(
        args:["sample_id"],
        kwargs:[
            "type": "test_sample",
            "barcode": null,
        ],
        name:"create_metamap",
    )
    return parser.parse_args(arguments)
}

/**
 * Take an input (file or directory) and return a channel of
 * named samples.
 *
 * @param input Top level input file or folder to locate sample data.
 * @param sample string to name single sample data.
 * @param sample_sheet Path to sample sheet CSV file.
 * @param sanitize regularize inputs from EPI2ME platform.
 * @param output output location, required if sanitize==true
 * @param min_barcode Minimum barcode to accept.
 * @param max_barcode Maximum (inclusive) barcode to accept.
 * @param input_type File input format. Default: fastq
 *
 * @return Channel of tuples (path, map(sample_id, type, barcode))
 */
def sample_ingress(Map arguments)
{
    def parser = new ArgumentParser(
        args:["input"],
        kwargs:[
            "sample":null, "sample_sheet":null, "sanitize":false, "output":null,
            "min_barcode":0, "max_barcode":Integer.MAX_VALUE, "input_type":"fastq"],
        name:"sample_ingress")
    Map margs = parser.parse_args(arguments)

    if (margs.sanitize && margs.output == null) {
        throw new Exception("Argument 'output' required if 'sanitize' is true.")
    }


    println("Checking input.")
    input = file(margs.input)

    // Handle file input
    if (input.isFile()) {
        // Assume sample is a string at this point
        println('Single file input detected.')
        if (margs.sample_sheet) {
            println('Warning: `--sample_sheet` given but single file input found. Ignoring.')
        }
        return handle_single_file(input, margs.sample)
    }

    // Handle directory input
    if (input.isDirectory()) {
        // EPI2ME harness
        if (margs.sanitize) {
            staging = file(margs.output).resolve("staging")
            input = sanitize_sample(input, staging, margs.input_type)
        }

        // Get barcoded and non barcoded subdirectories
        (barcoded, non_barcoded) = get_subdirectories(input)

        // Case 03: If no subdirectories, handle the single dir
        if (!barcoded && !non_barcoded) {
            println("Single directory input detected.")
            if (margs.sample_sheet) {
                println('Warning: `--sample_sheet` given but single non-barcode directory found. Ignoring.')
            }
            return handle_flat_dir(input, margs.sample, margs.input_type)
        }

        if (margs.sample) {
            println('Warning: `--sample` given but multiple directories found, ignoring.')
        }

        // Case 01, 02, 04: Handle barcoded and non_barcoded dirs
        // Handle barcoded folders
        barcoded_samples = Channel.empty()
        if (barcoded) {
            println("Barcoded directories detected.")
            sample_sheet = null
            if (margs.sample_sheet) {
                sample_sheet = get_sample_sheet(margs.sample_sheet)
            }
            barcoded_samples = handle_barcoded_dirs(barcoded, sample_sheet, margs.min_barcode, margs.max_barcode, margs.input_type)
        }

        non_barcoded_samples = Channel.empty()
        if (non_barcoded) {
            println("Non barcoded directories detected.")
            if (!barcoded && margs.sample_sheet) {
                println('Warning: `--sample_sheet` given but no barcode directories found.')
            }
            non_barcoded_samples = handle_non_barcoded_dirs(non_barcoded, margs.input_type)
        }

        return barcoded_samples.mix(non_barcoded_samples)
    }
}
