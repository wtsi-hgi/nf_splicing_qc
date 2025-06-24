/* ---- sequencing data QC pipeline ---- */

/* -- load modules -- */

/* -- load subworkflows -- */

/* -- define functions -- */

def helpMessage() {
    log.info """
Usage:
    nextflow run nf_splicing_qc/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet                path of the sample sheet
        --outdir                      the directory path of output results, default: the current directory
    
    Optional arguments:

    """
}

/* -- initialising parameters -- */
params.help                        = null
params.sample_sheet                = null
params.outdir                      = params.outdir                      ?: "$PWD"

/* -- check parameters -- */
if (params.help) {
    helpMessage()
    exit 0
}

// if (params.sample_sheet) {
//     ch_input = Channel.fromPath(file(params.sample_sheet), checkIfExists: true)
//                       .splitCsv(header: true, sep: ",")
//                       .map { row -> 
//                         def sample_id = "${row.sample}_${row.replicate}"
//                         tuple(sample_id, row.sample, row.replicate, row.directory, row.read1, row.read2, row.reference, row.barcode) }
// } else {
//     helpMessage()
//     log.info("Error: Please specify the full path of the sample sheet!\n")
//     exit 1
// }

if (!file(params.outdir).isDirectory()) {
    log.error("Invalid output directory: ${params.outdir}. Please specify a valid directory.")
    exit 1
}

/* -- check software exist -- */
// def required_tools = ['bwa', 'hisat2', 'samtools', 'bamtools', 'flash2', 'fastp']
// check_required(required_tools)


/* -- workflow -- */
workflow splicing_qc {
    
}
