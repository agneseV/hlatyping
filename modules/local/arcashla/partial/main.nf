process ARCASHLA_PARTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.6.0--hdfd78af_2':
        'biocontainers/arcas-hla:0.6.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads), path(genotype)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.partial_genotype.json"), emit: partial_genotype
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.partial_alignment.p"), optional: true, emit: alignment
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    def input = meta.single_end ? "--single ${reads[0]}" : "${reads[0]} ${reads[1]}"
    assert reads.size() == (meta.single_end ? 1 : 2) : 'Expected one FASTQ for SE or two for PE'
    """
    arcasHLA partial \
        ${args} \
        -G ${genotype} \
        -t ${task.cpus} \
        -o . \
        --temp . \
        ${input}

    shopt -s nullglob
    result=(*.partial_genotype.json)
    [[ \${#result[@]} -eq 1 ]] || { echo 'Expected exactly one partial genotype JSON' >&2; exit 1; }
    src="\${result[0]}"
    stem="\${src%.partial_genotype.json}"
    [[ "\$src" == '${prefix}.partial_genotype.json' ]] || mv "\$src" '${prefix}.partial_genotype.json'
    if [[ -f "\${stem}.partial_alignment.p" && "\${stem}" != '${prefix}' ]]; then
        mv "\${stem}.partial_alignment.p" '${prefix}.partial_alignment.p'
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcas-hla: '0.6.0'
    END_VERSIONS
    """
}
