#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { splicing_qc } from './workflows/splicing_qc.nf'

workflow {
    splicing_qc()
}
