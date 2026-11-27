//
// Align RNA reads: with STAR, then sort and mark duplicates (nf-core/oncoanalyser dev's route, rna_aligner 'star'), or
// with bwa-mem3 against the genome with transcript contigs appended, then lift alignments back to genome coordinates
// with TARS (rna_aligner 'bwa-mem3')
//

include { BWAMEM3_ALIGN_RNA    } from '../../../modules/local/bwa-mem3/mem/rna/main'
include { FASTP_SPLIT          } from '../../../modules/local/fastp/split/main'
include { GATK4_MARKDUPLICATES } from '../../../modules/local/gatk4/markduplicates/main'
include { SAMTOOLS_SORT        } from '../../../modules/nf-core/samtools/sort/main'
include { STAR_ALIGN           } from '../../../modules/local/star/align/main'
include { TARS                 } from '../../../modules/local/tars/main'

workflow READ_ALIGNMENT_RNA {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_fastq             // channel: [mandatory] [ meta, fastq_info, fastq_fwd, fastq_rev ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta (bwa-mem3)
    genome_version       // channel: [mandatory] genome version (bwa-mem3)
    genome_fai           // channel: [mandatory] /path/to/genome_fai (bwa-mem3)
    genome_dict          // channel: [mandatory] /path/to/genome_dict (bwa-mem3)
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/ (bwa-mem3)
    contigs_mapping_rna  // channel: [mandatory] /path/to/contigs_mapping_rna (bwa-mem3)
    unmap_regions_rna    // channel: [mandatory] /path/to/unmap_regions_rna (bwa-mem3)
    genome_star_index    // channel: [mandatory] /path/to/genome_star_index/ (star)

    // Params
    max_fastq_records    // numeric: [optional]  max number of FASTQ records per split (bwa-mem3)
    rna_aligner          //  string: [mandatory] 'star' or 'bwa-mem3'

    main:
    //
    // STEP: Handle inputs
    //
    // Sort inputs
    // runnable: channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    // skip: channel: [ meta ]

    // NOTE(LN): RNA is only ever sequenced for the tumor sample
    ch_inputs_rna_sorted = ch_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            def has_inputs = fastq_fwd && fastq_rev
            runnable: fastq_info.sample_type == 'tumor' && has_inputs
            skip: fastq_info.sample_type == 'tumor' && ! has_inputs
              return meta
        }

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_rna = channel.empty()
    // channel: [ meta, markdups_dir ]
    ch_outputs_rna_dir = channel.empty()

    if (rna_aligner == 'star') {

        // Create FASTQ input channel
        // NOTE(fork): STAR takes all lanes of a sample in one run, so FASTQs are not split here (max_fastq_records)
        // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
        ch_fastq_inputs = ch_inputs_rna_sorted.runnable
            .map { meta, fastq_info, fastq_fwd, fastq_rev ->

                // NOTE(SW): initial map sets defaults and conventional ordering of selected fields, merging then overwrites / adds while preserving order
                def rg_id = [fastq_info.sample_id, fastq_info.library_id, fastq_info.lane, fastq_info.flowcell].findAll().join('.')
                def rg_entries = [ID: rg_id, SM: fastq_info.sample_id, LB: fastq_info.library_id] + fastq_info.rg_fields
                def rg_line = rg_entries.collect { k, v -> "'${k}:${v}'" }.join(' ')

                def meta_fastq = [
                    key: meta.group_id,
                    id: "${meta.group_id}_${fastq_info.sample_id}",
                    sample_id: fastq_info.sample_id,
                    rg_line: rg_line,
                ]

                return [meta_fastq, fastq_fwd, fastq_rev]
            }

        //
        // MODULE: STAR alignment
        //
        // Create process input channel
        // First, count expected FASTQ pairs per sample for non-blocking groupTuple op
        // channel: [ meta_group, group_size ]
        ch_fastq_counts = ch_fastq_inputs
            .map { meta_fastq, _fastq_fwd, _fastq_rev ->
                def meta_group = [key: meta_fastq.key]
                return [meta_group, meta_fastq]
            }
            .groupTuple()
            .map { meta_group, meta_fastqs -> return [meta_group, meta_fastqs.size()] }

        // Now, group with expected size to proceed without blocking
        // channel: [ meta_star, [ rg_line, ... ], [ fastq_fwd, ... ], [ fastq_rev, ... ] ]
        ch_star_inputs = ch_fastq_counts
            // channel: [ [ meta_group, count ], [ meta_group, meta_fastq, fastq_fwd, fastq_rev ] ]
            .cross(
                // First element to match meta_group above for `cross`
                ch_fastq_inputs.map { meta_fastq, fastq_fwd, fastq_rev -> [[key: meta_fastq.key], meta_fastq, fastq_fwd, fastq_rev] }
            )
            // channel: [ GroupKey(meta_star, size), rg_line, fastq_fwd, fastq_rev ]
            .map { count_tuple, inputs_tuple ->
                def group_size = count_tuple[1]
                def (_meta_group, meta_fastq, fastq_fwd, fastq_rev) = inputs_tuple

                def meta_star = [
                    key: meta_fastq.key,
                    id: meta_fastq.id,
                    sample_id: meta_fastq.sample_id,
                ]

                return tuple(groupKey(meta_star, group_size), meta_fastq.rg_line, fastq_fwd, fastq_rev)
            }
            // channel: [ meta_star, rg_lines, fastq_fwds, fastq_revs ]
            .groupTuple()
            // Sort for stable inputs: lanes in read group order (the read group ID is unique per lane), not in the
            // order the FASTQ entries arrived
            .map { meta_star, rg_lines, fastq_fwds, fastq_revs ->
                def entries_sorted = [rg_lines, fastq_fwds, fastq_revs]
                    // [ [ rg_line, fastq_fwd, fastq_rev ], ... ]
                    .transpose()
                    .sort { a, b -> a[0] <=> b[0] ?: a[1].name <=> b[1].name }
                    // [ [ rg_line, ... ], [ fastq_fwd, ... ], [ fastq_rev, ... ] ]
                    .transpose()

                // Unpack cleanly
                def (rg_lines_sorted, fastq_fwds_sorted, fastq_revs_sorted) = entries_sorted

                return [meta_star, rg_lines_sorted, fastq_fwds_sorted, fastq_revs_sorted]
            }

        // Run process
        STAR_ALIGN(
            ch_star_inputs,
            genome_star_index,
        )

        //
        // MODULE: SAMtools sort
        //
        // Create process input channel
        // channel: [ meta_sort, aln ]
        ch_sort_inputs = channel.topic('star_align_bam')
            .map { meta_star, aln ->
                def meta_sort = meta_star + [prefix: meta_star.sample_id]
                return [meta_sort, aln]
            }

        // Run process
        SAMTOOLS_SORT(
            ch_sort_inputs,
        )

        //
        // MODULE: GATK4 MarkDuplicates
        //
        // Create process input channel, restoring original meta to source the RNA sample name
        // channel: [ meta_markdups, aln ]
        ch_markdups_inputs = WorkflowOncoanalyser.restoreMeta(channel.topic('samtools_sort_bam'), ch_inputs)
            .map { meta, aln ->

                def sample_id = Utils.getTumorRnaSampleName(meta)

                def meta_markdups = [
                    key: meta.group_id,
                    id: "${meta.group_id}_${sample_id}",
                    sample_id: sample_id,
                ]

                return [meta_markdups, aln]
            }

        // Run process
        GATK4_MARKDUPLICATES(
            ch_markdups_inputs,
        )

        //
        // STEP: Handle outputs
        //
        // NOTE(fork): the MarkDuplicates directory (<sample_id>.md.bam) takes the place of the RNA REDUX directory in
        // the calling workflow, so no alignments are passed on to REDUX here
        // channel: [ meta, [], [] ]
        ch_outputs_rna = ch_inputs.map { meta -> [meta, [], []] }

        // Set outputs, restoring original meta
        // channel: [ meta, markdups_dir ]
        ch_outputs_rna_dir = channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(channel.topic('gatk4_markduplicates_dir'), ch_inputs),
                ch_inputs_rna_sorted.skip.unique().map { meta -> [meta, []] },
            )

    } else {

        // Create FASTQ input channel
        // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
        ch_fastq_inputs = channel.empty()
            .mix(ch_inputs_rna_sorted.runnable)
            .map { meta, fastq_info, fastq_fwd, fastq_rev ->

                // NOTE(SW): initial map sets defaults and conventional ordering of selected fields, merging then overwrites / adds while preserving order
                def rg_id = [fastq_info.sample_id, fastq_info.library_id, fastq_info.lane, fastq_info.flowcell].findAll().join('.')
                def rg_entries = [ID: rg_id, SM: fastq_info.sample_id, LB: fastq_info.library_id] + fastq_info.rg_fields
                def rg_line = '@RG\\t' + rg_entries.collect { k, v -> "${k}:${v}" }.join('\\t')

                def meta_fastq = [
                    key: meta.group_id,
                    id: "${meta.group_id}_${fastq_info.sample_id}_${fastq_info.library_id}_${fastq_info.lane}",
                    rg_line: rg_line,
                    sample_id: fastq_info.sample_id,
                    library_id: fastq_info.library_id,
                    lane: fastq_info.lane,
                    output_file_id: rg_id,
                    sample_type: fastq_info.sample_type,
                ]

                if (fastq_info.flowcell) {
                    meta_fastq.id = "${meta_fastq.id}_${fastq_info.flowcell}"
                }

                return [meta_fastq, fastq_fwd, fastq_rev]

            }

        //
        // MODULE: fastp
        //
        // Split FASTQ into chunks if requested for distributed processing
        // channel: [ meta_fastq_ready, fastq_fwd, fastq_fwd ]
        ch_fastqs_ready = channel.empty()
        // NOTE(SW): required for strict syntax without params block declaration
        if (max_fastq_records.toInteger() > 0) {

            // Run process
            FASTP_SPLIT(
                ch_fastq_inputs,
                // NOTE(SW): required for strict syntax without params block declaration
                max_fastq_records.toInteger(),
            )

            // NOTE(LN): the transpose operator pairs the R1 and R2 chunks by index, and also covers the single chunk case
            // where fastp emits one file per read rather than a list
            ch_fastqs_ready = channel.topic('fastp_split_fastq')
                .transpose()
                .map { meta_fastq, fwd, rev ->

                    def split_fwd = fwd.name.replaceAll('\\..+$', '')
                    def split_rev = rev.name.replaceAll('\\..+$', '')

                    assert split_fwd == split_rev

                    // NOTE(SW): split allows meta_fastq_ready to be unique, which is required during reunite below
                    def meta_fastq_ready = meta_fastq + [id: "${meta_fastq.id}_${split_fwd}", split: split_fwd]

                    return [meta_fastq_ready, fwd, rev]
                }

        } else {

            ch_fastqs_ready = ch_fastq_inputs
                .map { meta_fastq, fastq_fwd, fastq_rev ->

                    def meta_fastq_ready = meta_fastq + [split: null]

                    return [meta_fastq_ready, fastq_fwd, fastq_rev]
                }

        }

        //
        // MODULE: BWA-MEM3
        //
        // Create process input channel
        // channel: [ meta_bwamem3, fastq_fwd, fastq_rev ]
        ch_bwamem3_inputs = ch_fastqs_ready
            .map { meta_fastq_ready, fastq_fwd, fastq_rev ->
                def meta_bwamem3 = meta_fastq_ready.clone()
                return [meta_bwamem3, fastq_fwd, fastq_rev]
            }

        // Run process
        BWAMEM3_ALIGN_RNA(
            ch_bwamem3_inputs,
            genome_fasta,
            genome_bwamem2_index,
        )

        // Reunite BAMs
        // Count expected BAMs per sample for non-blocking groupTuple op
        // channel: [ meta_group, group_size ]
        ch_sample_fastq_counts = ch_bwamem3_inputs
            .map { meta_bwamem3, _reads_fwd, _reads_rev ->

                def meta_group = [
                    key: meta_bwamem3.key,
                    sample_type: meta_bwamem3.sample_type,
                ]

                return [meta_group, meta_bwamem3]
            }
            .groupTuple()
            .map { meta_group, metas_bwamem3 -> return [meta_group, metas_bwamem3.size()] }

        // Now, group with expected size
        // NOTE(LN): RNA alignments are name-grouped and therefore unindexed, so no index is carried here
        // channel: [ meta_group, [aln, ...] ]
        ch_alns_united = ch_sample_fastq_counts
            // channel: [ [ meta_group, count ], [ meta_group, aln ] ]
            .cross(
                // First element to match meta_group above for `cross`
                channel.topic('bwamem3_align_rna_bam').map { meta_bwamem3, aln -> [[key: meta_bwamem3.key, sample_type: meta_bwamem3.sample_type], aln] }
            )
            .map { count_tuple, inputs_tuple ->
                def group_size = count_tuple[1]
                def (meta_group, aln) = inputs_tuple

                return tuple(groupKey(meta_group, group_size), aln)
            }
            // Lane BAMs in file-name order, not in the order the lanes finished (see read_alignment_dna)
            .groupTuple(sort: { a, b -> a.name <=> b.name })

        //
        // MODULE: TARS
        //
        // Create process input channel, restoring original meta to source the RNA sample name
        // channel: [ meta_tars, [aln, ...] ]
        ch_tars_inputs = WorkflowOncoanalyser.restoreMeta(ch_alns_united, ch_inputs)
            .map { meta, alns ->

                def sample_id = Utils.getTumorRnaSampleName(meta)

                def meta_tars = [
                    key: meta.group_id,
                    id: "${meta.group_id}_${sample_id}",
                    sample_id: sample_id,
                ]

                return [meta_tars, alns]
            }

        // Run process
        TARS(
            ch_tars_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            genome_dict,
            contigs_mapping_rna,
            unmap_regions_rna,
        )

        //
        // STEP: Handle outputs
        //
        // Set outputs, restoring original meta
        // NOTE(LN): Tars emits one BAM per sample, but it is carried as a list so that the output shape matches the DNA
        // read alignment subworkflow
        // channel: [ meta, [aln, ...], [idx, ...] ]
        ch_outputs_rna = channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(channel.topic('tars_bam'), ch_inputs)
                    .map { meta, aln, idx -> [meta, [aln], [idx]] },
                ch_inputs_rna_sorted.skip.unique().map { meta -> [meta, [], []] },
            )

        // No MarkDuplicates directory: the RNA REDUX directory is used
        // channel: [ meta, [] ]
        ch_outputs_rna_dir = ch_inputs.map { meta -> [meta, []] }

    }

    emit:
    rna     = ch_outputs_rna     // channel: [ meta, [aln, ...], [idx, ...] ]
    rna_dir = ch_outputs_rna_dir // channel: [ meta, markdups_dir ]
}
