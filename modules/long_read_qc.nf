process filtlong {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_RL.filtered.fastq.gz"),    emit: filtered_reads
    tuple val(sample_id), path("${sample_id}_filtlong_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: filtlong\\n"                                          >> ${sample_id}_filtlong_provenance.yml
    printf -- "  tools:\\n"                                                          >> ${sample_id}_filtlong_provenance.yml
    printf -- "    - tool_name: filtlong\\n"                                         >> ${sample_id}_filtlong_provenance.yml
    printf -- "      tool_version: \$(filtlong --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_filtlong_provenance.yml
    printf -- "      parameters:\\n"                                                 >> ${sample_id}_filtlong_provenance.yml
    printf -- "        - parameter: --min_length\\n"                                 >> ${sample_id}_filtlong_provenance.yml
    printf -- "          value: ${params.filtlong_min_length}\\n"                    >> ${sample_id}_filtlong_provenance.yml
    printf -- "        - parameter: --keep_percent\\n"                               >> ${sample_id}_filtlong_provenance.yml
    printf -- "          value: ${params.filtlong_keep_percent}\\n"                  >> ${sample_id}_filtlong_provenance.yml

    filtlong \
      --min_length   ${params.filtlong_min_length} \
      --keep_percent ${params.filtlong_keep_percent} \
      ${reads} | \
        gzip > ${sample_id}_RL.filtered.fastq.gz
    """
}

process nanoq {

    tag { sample_id + ' / ' + pre_or_post_filter }

    input:
    tuple val(sample_id), path(reads), val(pre_or_post_filter)

    output:
    tuple val(sample_id), path("${sample_id}_nanoq_*.csv"),                                emit: report
    tuple val(sample_id), path("${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: nanoq_${pre_or_post_filter}\\n"                    >> ${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml
    printf -- "    - tool_name: nanoq\\n"                                         >> ${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml
    printf -- "      tool_version: \$(nanoq --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml
    printf -- "      parameters:\\n"                                              >> ${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml
    printf -- "        - parameter: --stats\\n"                                   >> ${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml
    printf -- "          value: null\\n"                                          >> ${sample_id}_nanoq_${pre_or_post_filter}_provenance.yml

    nanoq --header --stats --input ${reads} | tr ' ' ',' > ${sample_id}_nanoq_${pre_or_post_filter}.csv
    """
}

process merge_nanoq_reports {

    tag { sample_id }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_nanoq.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(nanoq_pre_filter), path(nanoq_post_filter)

    output:
    tuple val(sample_id), path("${sample_id}_nanoq.csv")

    script:
    """
    merge_nanoq_reports.py --sample-id ${sample_id} --pre-filter ${nanoq_pre_filter} --post-filter ${nanoq_post_filter} > ${sample_id}_nanoq.csv
    """
}


process bandage {

    tag { sample_id + ' / ' + assembly_mode }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_${assembler}_${assembly_mode}_bandage.png", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly_graph), val(assembler), val(assembly_mode)

    output:
    tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_bandage.png")

    script:
    """
    Bandage image ${assembly_graph} ${sample_id}_${assembler}_${assembly_mode}_bandage.png
    """
}
