process NEO_ANNOTATE_FUSIONS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-isofox:2.0.1--hdfd78af_0' :
        'biocontainers/hmftools-isofox:2.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(neo_finder_dir), path(aln), path(idx)
    val read_length
    path genome_fasta
    val genome_ver
    path genome_fai
    path ensembl_data_resources

    output:
    tuple val(meta), path('*isf.neoepitope.tsv')                    , topic: neo_annotated_fusions_tsv
    tuple val(meta), val('neo_annotate_fusions'), path('.command.*'), topic: command_files
    path 'versions.yml'                                             , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    // ISOFOX 2.0 (the STAR route, --rna_aligner star) has no -tumor and reads NEO's file by -sample, the DNA sample, as
    // upstream nf-core/oncoanalyser runs it; NEO's scorer reads the result by the RNA sample, so it is renamed after
    def star_route = params.rna_aligner == 'star'
    def sample_args = star_route ? "-sample ${meta.tumor_dna_id}" : "-sample ${meta.tumor_rna_id} -tumor ${meta.tumor_dna_id}"
    def rename_cmd = star_route && meta.tumor_dna_id != meta.tumor_rna_id ? "mv ${meta.tumor_dna_id}.isf.neoepitope.tsv ${meta.tumor_rna_id}.isf.neoepitope.tsv" : ''

    """
    isofox \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        ${sample_args} \\
        -bam_file ${aln} \\
        -functions NEO_EPITOPES \\
        -read_length ${read_length} \\
        -neo_dir ${neo_finder_dir} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir ./

    ${rename_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox ${args} -version | sed -n '/^Isofox version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_rna_id}.isf.neoepitope.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
