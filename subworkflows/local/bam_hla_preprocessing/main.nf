//
// Convert a genome-aligned BAM into HLA-enriched FASTQ reads.
//
// The workflow keeps primary, non-duplicate reads overlapping the supplied HLA BED
// together with their mates, adds fully-unmapped read pairs, merges both BAM subsets,
// and converts the result to non-interleaved FASTQ files.
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHECK_PAIRED                              } from '../../../modules/local/check_paired'
include { SAMTOOLS_INDEX                            } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HLA        } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_UNMAPPED   } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_MERGE                            } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_COLLATEFASTQ                     } from '../../../modules/nf-core/samtools/collatefastq/main'
include { SAMTOOLS_BEDCOV                           } from '../../../modules/local/samtools/bedcov/main'
include { COV_CHECK                                 } from '../../../modules/local/cov_check/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BAM_HLA_PREPROCESSING {
    take:
    ch_bam // channel: [ val(meta), path(bam) ]
    ch_bed // value channel: [ val(meta), path(bed) ]

    main:
    // The main workflow carries BAMs as a one-element file list; normalize to a single Path.
    ch_bam
        .map { meta, files -> [meta, files instanceof List ? files[0] : files] }
        .set { ch_bam_input }

    // Detect an existing index alongside each original input BAM
    ch_existing_bai = ch_bam_input.map { meta, bam ->

        def candidates = [
            file("${bam}.bai"),
            file(bam.toString().replaceFirst(/\.bam$/, '.bai'))
        ]

        def bai = candidates.find { it.exists() }

        tuple(meta.id, bai ?: [])
    }


    // Determine endedness because SAMTOOLS_COLLATEFASTQ uses meta.single_end.
    CHECK_PAIRED(ch_bam_input)
    CHECK_PAIRED.out.reads
        .map { meta, bam, single_end ->
            [meta + [single_end: single_end.text.toBoolean()], bam]
        }
        .set { ch_bam_pe_corrected }

    // Match the corrected BAM metadata to the original index
    ch_bam_with_bai = ch_bam_pe_corrected
        .map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(ch_existing_bai)
        .map { id, meta, bam, bai -> tuple(meta, bam, bai) }

    // Separate indexed and non-indexed BAMs
    ch_index_branch = ch_bam_with_bai.branch { meta, bam, bai ->
        indexed: bai != []
        needs_index: true
    }

    // Index only BAMs without an existing index
    ch_to_index = ch_index_branch.needs_index
        .map { meta, bam, bai -> tuple(meta, bam) }

    SAMTOOLS_INDEX(ch_to_index)

    // Combine existing and newly generated indexes
    ch_bam_indexed = ch_index_branch.indexed
        .mix(
            ch_to_index.join(SAMTOOLS_INDEX.out.index)
        )

    // HLA/MHC-overlapping primary, non-duplicate reads plus their mates.
    // -L is supplied by the SAMTOOLS_VIEW bed input; filtering flags are set in modules.config.
    SAMTOOLS_VIEW_HLA(
        ch_bam_indexed,
        [[:], [], []],
        [[:], []],
        ch_bed,
        '',
    )

    // Fully-unmapped primary, non-duplicate read pairs.
    SAMTOOLS_VIEW_UNMAPPED(
        ch_bam_indexed,
        [[:], [], []],
        [[:], []],
        [[:], []],
        '',
    )

    // Both BAMs derive from the same source BAM and are disjoint by construction.
    SAMTOOLS_VIEW_HLA.out.bam
        .join(SAMTOOLS_VIEW_UNMAPPED.out.bam)
        .map { meta, mhc_bam, unmapped_bam ->
            [meta, [mhc_bam, unmapped_bam], []]
        }
        .set { ch_bam_to_merge }

    SAMTOOLS_MERGE(
        ch_bam_to_merge,
        [[:], [], [], []],
    )

    // Match the existing hlatyping BAM -> FASTQ behaviour: separate R1/R2 files.
    def interleave = false
    SAMTOOLS_COLLATEFASTQ(
        SAMTOOLS_MERGE.out.bam,
        SAMTOOLS_MERGE.out.bam.map { meta, _bam -> [[id: meta.id], [], []] },
        interleave,
    )


    // Exon coverage test

    ch_hlabed = Channel.value(
        file(params.hla_exon_bed, checkIfExists: true)
    )

    SAMTOOLS_BEDCOV(
        ch_bam_indexed,
        ch_hlabed
    )

    COV_CHECK(
    SAMTOOLS_BEDCOV.out.coverage
    )

    emit:
    reads         = SAMTOOLS_COLLATEFASTQ.out.fastq
    singleton     = SAMTOOLS_COLLATEFASTQ.out.fastq_singleton
    other         = SAMTOOLS_COLLATEFASTQ.out.fastq_other
    merged_bam    = SAMTOOLS_MERGE.out.bam
    mhc_bam       = SAMTOOLS_VIEW_HLA.out.bam
    unmapped_bam  = SAMTOOLS_VIEW_UNMAPPED.out.bam
    cov_check     = COV_CHECK.out.note    

}
