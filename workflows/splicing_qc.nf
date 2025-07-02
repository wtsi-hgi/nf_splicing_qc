/* ---- sequencing data QC pipeline ---- */

/* -- load modules -- */

/* -- load subworkflows -- */

/* -- define functions -- */

def helpMessage() {
    log.info """
Usage:
    nextflow run nf_splicing_qc/main.nf --sample_sheet "/path/of/sample/sheet" --experiment_design "/path/of/experiment/design"

    Mandatory arguments:
    Basic:
        --sample_sheet                path of the sample sheet
        --outdir                      the directory path of output results, default: the current directory
    
    DiMSum:
        --experiment_design           path of the experiment design file

    Optional arguments:
    Cutadapt:
        --minimum-length               discard reads shorter than this length after trimming, default: 38
        --error-rate                   maximum allowed error rate as value between 0 and 1, default: 0.2
        --overlap                      require minimal overlap between read and adapter for an adapter to be found, default: 3

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

if (params.sample_sheet) {
    def sep = params.sample_sheet.endsWith('.tsv') ? '\t' : ','
    ch_input = Channel.fromPath(file(params.sample_sheet), checkIfExists: true)
                      .splitCsv(header: true, sep: sep)
                      .map { row -> 
                        def sample_id = "${row.sample}_${row.replicate}"
                        tuple(sample_id, row.sample, row.replicate, row.directory, row.read1, row.read2, cutadapt_g, cutadapt_G) }
} else {
    helpMessage()
    log.info("Error: Please specify the full path of the sample sheet!\n")
    exit 1
}

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
