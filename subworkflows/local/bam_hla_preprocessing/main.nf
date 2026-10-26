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

    // Determine endedness because SAMTOOLS_COLLATEFASTQ uses meta.single_end.
    CHECK_PAIRED(ch_bam_input)
    CHECK_PAIRED.out.reads
        .map { meta, bam, single_end ->
            [meta + [single_end: single_end.text.toBoolean()], bam]
        }
        .set { ch_bam_pe_corrected }

    // --fetch-pairs (-P) requires an indexed regular BAM file.
    SAMTOOLS_INDEX(ch_bam_pe_corrected)
    ch_bam_pe_corrected
        .join(SAMTOOLS_INDEX.out.index)
        .set { ch_bam_indexed }

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

    emit:
    reads       = SAMTOOLS_COLLATEFASTQ.out.fastq
    singleton   = SAMTOOLS_COLLATEFASTQ.out.fastq_singleton
    other       = SAMTOOLS_COLLATEFASTQ.out.fastq_other
    merged_bam  = SAMTOOLS_MERGE.out.bam
    mhc_bam     = SAMTOOLS_VIEW_HLA.out.bam
    unmapped_bam = SAMTOOLS_VIEW_UNMAPPED.out.bam
}
