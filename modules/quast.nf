process quast {

    tag { sample_id + ' / ' + assembly_mode }

    input:
      tuple val(sample_id), path(assembly), val(assembler), val(assembly_mode)

    output:
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_quast.tsv"), val(assembler), val(assembly_mode), emit: tsv
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml"),                          emit: provenance

    script:
      """
      printf -- "- process_name: quast\\n"                                                 >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "  tools:\\n"                                                              >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "    - tool_name: quast\\n"                                                >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "      tool_version: \$(quast --version | cut -d ' ' -f 2 | tr -d 'v')\\n" >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "      parameters:\\n"                                                     >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "        - parameter: --space-efficient\\n"                                >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "          value: null\\n"                                                 >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "        - parameter: --fast\\n"                                           >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml
      printf -- "          value: null\\n"                                                 >> ${sample_id}_${assembler}_${assembly_mode}_quast_provenance.yml

      quast \
        --threads ${task.cpus} \
        --space-efficient \
        --fast \
        --output-dir ${sample_id} \
        ${assembly}

      mv ${sample_id}/transposed_report.tsv ${sample_id}_${assembler}_${assembly_mode}_quast.tsv
      """
}

process parse_quast_report {

    tag { sample_id + ' / ' + assembly_mode }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_${assembler}_${assembly_mode}_quast.csv", mode: 'copy'

    input:
      tuple val(sample_id), path(quast_report), val(assembler), val(assembly_mode)

    output:
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_quast.csv")

    script:
      """
      parse_quast_report.py ${quast_report} > ${sample_id}_${assembler}_${assembly_mode}_quast.csv
      """
}
