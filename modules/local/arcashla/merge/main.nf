
process ARCASHLA_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.6.0--hdfd78af_2' :
        'biocontainers/arcas-hla:0.6.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(genotypes)

    output:
    tuple val(meta), path("*.genotypes.tsv"), emit: genotypes
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    mkdir -p input

    cp *.genotype.json input/

    arcasHLA merge \\
        ${args} \\
        --run ${prefix} \\
        --indir input \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcasHLA: 0.6.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.genotypes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcasHLA: 0.6.0
    END_VERSIONS
    """
}
