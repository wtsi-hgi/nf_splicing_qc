workflow check_input_files {
    take:
    ch_sample

    main:


    emit:

}

process check_files {
    label 'process_single'

    input:

    output:

    script:
}

/*
sample_sheet:
sample, rep, read_dir, read1, read2

experiment_design:
sample_name, experiment_replicate, selection_id, selection_replicate, technical_replicate, pair1, pair2, cutadapt5First, cutadapt5Second
*/