process COV_CHECK {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    tuple val(meta), path(exon_coverage)

    output:
    tuple val(meta), path('*.cov_check.txt'), emit: note

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hla_cov_check.py \\
        ${exon_coverage} \\
        --output ${prefix}.cov_check.txt \\
        --summary ${prefix}.gene_coverage.tsv \\
        ${args}
    """
}

