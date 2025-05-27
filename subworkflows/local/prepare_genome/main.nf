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
    gens_interval_list = Channel.empty()
    mutect2_target_bed = Channel.empty()
    versions = Channel.empty()

    // Only run GATK4_CREATESEQUENCEDICTIONARY and generate dict if no user_dict is provided
    fasta_for_dict = fasta
        .join(user_dict, remainder: true)
        .filter { _meta, _fasta, dict_ -> !dict_ }
        .map { meta, fasta_, _dict -> [meta, fasta_] }

    GATK4_CREATESEQUENCEDICTIONARY(fasta_for_dict)

    dict = user_dict.mix(GATK4_CREATESEQUENCEDICTIONARY.out.dict).collect()

    // Only run SAMTOOLS_FAIDX and generate fai if no user_fai is provided
    fasta_for_fai = fasta
        .join(user_fai, remainder: true)
        .filter { _meta, _fasta, fai_ -> !fai_ }
        .map { meta, fasta_, _fai -> [meta, fasta_] }

    SAMTOOLS_FAIDX(fasta_for_fai, [[:], []])

    fai = user_fai.mix(SAMTOOLS_FAIDX.out.fai).collect()

    // Only run GATK4_PREPROCESSINTERVALS_GENS and generate gens_interval_list if no user_gens_interval_list is provided
    fasta_for_interval_list = fasta
        .join(user_gens_interval_list, remainder: true)
        .filter { _meta, _fasta, interval_list_ -> !interval_list_ }
        .map { meta, fasta_, _interval_list -> [meta, fasta_] }

    GATK4_PREPROCESSINTERVALS_GENS(fasta_for_interval_list, fai.collect(), dict.collect(), [[:], []], [[:], []])

    gens_interval_list = user_gens_interval_list.mix(GATK4_PREPROCESSINTERVALS_GENS.out.interval_list).collect()

    // Only run BUILD_INTERVALS and generate mutect2_target_bed if no user_mutect2_target_bed is provided
    fai_for_intervals = fai
        .join(user_mutect2_target_bed, remainder: true)
        .filter { _meta, _fai, mutect2_target_bed_ -> !mutect2_target_bed_ }
        .map { meta, fai_, _mutect2_target_bed -> [meta, fai_] }

    BUILD_INTERVALS(fai_for_intervals, [], false)

    mutect2_target_bed = user_mutect2_target_bed.mix(BUILD_INTERVALS.out.output).collect()

    versions = versions.mix(BUILD_INTERVALS.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(GATK4_PREPROCESSINTERVALS_GENS.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    dict               // channel: [mandatory] [ val(meta), path(dict) ]
    fai                // channel: [mandatory] [ val(meta), path(fai) ]
    gens_interval_list // channel: [mandatory] [ val(meta), path(gens_interval_list) ]
    mutect2_target_bed // channel: [mandatory] [ val(meta), path(mutect2_target_bed) ]
    versions           // channel: path(versions.yml)
}
