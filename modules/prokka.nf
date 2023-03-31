process prokka {

    tag { sample_id + ' / ' + assembly_mode}

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_${assembler}_prokka.{gbk,gff}", mode: 'copy'

    input:
      tuple val(sample_id), path(assembly), val(assembler), val(assembly_mode)

    output:
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_prokka.gbk"),            emit: gbk
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_prokka.gff"),            emit: gff
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml"), emit: provenance

    script:
      """
      printf -- "- process_name: prokka\\n"                                          >> ${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml
      printf -- "  tools:\\n"                                                        >> ${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml
      printf -- "    - tool_name: prokka\\n"                                         >> ${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml  
      printf -- "      tool_version: \$(prokka --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml  
      printf -- "      parameters:\\n"                                               >> ${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml
      printf -- "        - parameter: --compliant\\n"                                >> ${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml
      printf -- "          value: null\\n"                                           >> ${sample_id}_${assembler}_${assembly_mode}_prokka_provenance.yml

      prokka --cpus ${task.cpus} --compliant --locustag ${sample_id} --centre "BCCDC-PHL" --prefix "${sample_id}" ${assembly}

      cp ${sample_id}/${sample_id}.gbk ${sample_id}_${assembler}_${assembly_mode}_prokka.gbk
      cp ${sample_id}/${sample_id}.gff ${sample_id}_${assembler}_${assembly_mode}_prokka.gff
      """
}
