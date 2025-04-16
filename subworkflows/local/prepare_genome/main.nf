include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_PREPROCESSINTERVALS      } from '../../../modules/nf-core/gatk4/preprocessintervals'

//  Prepare references
workflow PREPARE_GENOME {
    take:
    fasta                   // channel: [mandatory] [ val(meta), path(fasta) ]
    user_dict               // channel: [optional]  [ val(meta), path(dict) ]
    user_fai                // channel: [optional]  [ val(meta), path(fai) ]
    user_gens_interval_list // channel: [optional]  [ val(meta), path(interval_list) ]
    tools                   //   array: [mandatory] [ tools ]

    main:
    dict = Channel.empty()
    fai = Channel.empty()
    interval_list = Channel.empty()
    versions = Channel.empty()

    // If more than one file, then it means that the user has provided a fai file
    // So we can pass out a null channel and SAMTOOLS_FAIDX won't be run
    fasta_for_fai = fasta
        .mix(user_fai)
        .groupTuple()
        .map { meta, files ->
            files[1] ? null : [meta, files[0]]
        }

    // If more than one file, then it means that the user has provided a dict file
    // So we can pass out a null channel and GATK4_CREATESEQUENCEDICTIONARY won't be run
    fasta_for_dict = fasta
        .mix(user_dict)
        .groupTuple()
        .map { meta, files ->
            files[1] ? null : [meta, files[0]]
        }

    SAMTOOLS_FAIDX(fasta_for_fai, [[:], []])
    GATK4_CREATESEQUENCEDICTIONARY(fasta_for_dict)

    dict = user_dict.mix(GATK4_CREATESEQUENCEDICTIONARY.out.dict).collect()

    fai = user_fai.mix(SAMTOOLS_FAIDX.out.fai).collect()


    // If more than one file, then it means that the user has provided an interval list file
    // So we can pass out a null channel and GATK4_PREPROCESSINTERVALS won't be run

    fasta_for_interval_list = fasta
        .mix(user_gens_interval_list)
        .groupTuple()
        .map { meta, files ->
            files[1] || !tools.split(',').contains('gens') ? null : [meta, files[0]]
        }

    GATK4_PREPROCESSINTERVALS(fasta_for_interval_list, fai.collect(), dict.collect(), [[:], []], [[:], []])

    interval_list = user_gens_interval_list.mix(GATK4_PREPROCESSINTERVALS.out.interval_list).collect()

    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    dict
    fai
    interval_list
    versions
}
