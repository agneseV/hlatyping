process ARCASHLA_GENOTYPE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.6.0--hdfd78af_2':
        'biocontainers/arcas-hla:0.6.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.genotype.json"), emit: genotype
    tuple val(meta), path("*.em.json"), optional: true, emit: em
    tuple val(meta), path("*.alignment.p"), optional: true, emit: alignment
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    def input = meta.single_end ? "--single ${reads[0]}" : "${reads[0]} ${reads[1]}"
    assert reads.size() == (meta.single_end ? 1 : 2) : 'Expected one FASTQ for SE or two for PE'
    """
    arcasHLA genotype \
        ${args} \
        -t ${task.cpus} \
        -o . \
        --temp . \
        ${input}

    shopt -s nullglob
    result=(*.genotype.json)
    [[ \${#result[@]} -eq 1 ]] || { echo 'Expected exactly one genotype JSON' >&2; exit 1; }
    src="\${result[0]}"
    stem="\${src%.genotype.json}"
    [[ "\$src" == '${prefix}.genotype.json' ]] || mv "\$src" '${prefix}.genotype.json'
    for suffix in em.json alignment.p; do
        if [[ -f "\${stem}.\${suffix}" && "\${stem}.\${suffix}" != '${prefix}'".\${suffix}" ]]; then
            mv "\${stem}.\${suffix}" '${prefix}'".\${suffix}"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcas-hla: '0.6.0'
    END_VERSIONS
    """
}
