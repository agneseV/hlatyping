process SAMTOOLS_BEDCOV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path('*.exon_coverage.tsv'), emit: coverage
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-j'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Validate exactly four BED columns (chrom, start, end, name).
    awk 'BEGIN { FS=OFS="\\t" }
         /^#/ || /^track/ || /^browser/ { next }
         NF != 4 || \$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ || \$3 <= \$2 {
             print "ERROR: expected BED4 with nonempty, positive-length intervals at line " NR > "/dev/stderr"
             exit 1
         }' ${bed}

    # Use the explicitly staged index. -X accepts the BAM and index paths.
    samtools bedcov ${args} -X ${bed} ${bam} ${bai} > bedcov_raw.tsv

    awk 'BEGIN { FS=OFS="\\t"; print "chrom", "start", "end", "region", "depth_sum", "mean_depth" }
         NF != 5 { print "ERROR: unexpected samtools bedcov output" > "/dev/stderr"; exit 1 }
         { len=\$3-\$2; if (len <= 0) { print "ERROR: invalid BED interval" > "/dev/stderr"; exit 1 }
           printf "%s\\t%s\\t%s\\t%s\\t%s\\t%.6f\\n", \$1,\$2,\$3,\$4,\$5,\$5/len }
        ' bedcov_raw.tsv > ${prefix}.exon_coverage.tsv
    """
}
