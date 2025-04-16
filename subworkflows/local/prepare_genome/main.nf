include { GATK4_CREATESEQUENCEDICTIONARY                              } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_PREPROCESSINTERVALS as GATK4_PREPROCESSINTERVALS_GENS } from '../../../modules/nf-core/gatk4/preprocessintervals'
include { GAWK as BUILD_INTERVALS                                     } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX                                              } from '../../../modules/nf-core/samtools/faidx'

//  Prepare references
workflow PREPARE_GENOME {
    take:
    fasta                   // channel: [mandatory] [ val(meta), path(fasta) ]
    user_dict               // channel: [optional]  [ val(meta), path(dict) ]
    user_fai                // channel: [optional]  [ val(meta), path(fai) ]
    user_gens_interval_list // channel: [optional]  [ val(meta), path(gens_interval_list) ]
    user_mutect2_target_bed // channel: [optional]  [ val(meta), path(mutect2_target_bed) ]
    tools                   //   array: [mandatory] [ tools ]

    main:
    dict = Channel.empty()
    fai = Channel.empty()
    interval_list = Channel.empty()
    versions = Channel.empty()

    // If more than one file, then it means that the user has provided a dict file
    // So we can pass out a null channel and GATK4_CREATESEQUENCEDICTIONARY won't be run
    fasta_for_dict = fasta
        .mix(user_dict)
        .groupTuple()
        .map { meta, files ->
            files[1] ? null : [meta, files[0]]
        }

    GATK4_CREATESEQUENCEDICTIONARY(fasta_for_dict)

    dict = user_dict.mix(GATK4_CREATESEQUENCEDICTIONARY.out.dict).collect()

    // If more than one file, then it means that the user has provided a fai file
    // So we can pass out a null channel and SAMTOOLS_FAIDX won't be run
    fasta_for_fai = fasta
        .mix(user_fai)
        .groupTuple()
        .map { meta, files ->
            files[1] ? null : [meta, files[0]]
        }

    SAMTOOLS_FAIDX(fasta_for_fai, [[:], []])

    fai = user_fai.mix(SAMTOOLS_FAIDX.out.fai).collect()

    // If more than one file, then it means that the user has provided an interval list file
    // So we can pass out a null channel and GATK4_PREPROCESSINTERVALS_GENS won't be run

    fasta_for_interval_list = fasta
        .mix(user_gens_interval_list)
        .groupTuple()
        .map { meta, files ->
            files[1] || !tools.split(',').contains('gens') ? null : [meta, files[0]]
        }

    GATK4_PREPROCESSINTERVALS_GENS(fasta_for_interval_list, fai.collect(), dict.collect(), [[:], []], [[:], []])

    gens_interval_list = user_gens_interval_list.mix(GATK4_PREPROCESSINTERVALS_GENS.out.interval_list).collect()

    // If more than one file, then it means that the user has provided a fai file
    // So we can pass out a null channel and SAMTOOLS_FAIDX won't be run
    fai_for_intervals = fai
        .mix(user_mutect2_target_bed)
        .groupTuple()
        .map { meta, files ->
            files[1] ? null : [meta, files[0]]
        }

    BUILD_INTERVALS(fai_for_intervals, [], false)

    mutect2_target_bed = user_mutect2_target_bed.mix(BUILD_INTERVALS.out.output).collect()

    versions = versions.mix(BUILD_INTERVALS.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(GATK4_PREPROCESSINTERVALS_GENS.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    dict
    fai
    gens_interval_list
    mutect2_target_bed
    versions
}
