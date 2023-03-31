process dragonflye {

    tag { sample_id + ' / ' + assembly_mode }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_dragonflye_${assembly_mode}.{fa,gfa,log}", mode: 'copy'

    input:
      tuple val(sample_id), path(reads), val(assembly_mode)

    output:
      tuple val(sample_id), path("${sample_id}_dragonflye_${assembly_mode}.fa"),  val("dragonflye"), val(assembly_mode), emit: assembly
      tuple val(sample_id), path("${sample_id}_dragonflye_${assembly_mode}_unpolished.gfa"), val("dragonflye"), val(assembly_mode), emit: assembly_graph
      tuple val(sample_id), path("${sample_id}_dragonflye_${assembly_mode}.log"), val("dragonflye"), val(assembly_mode), emit: log
      tuple val(sample_id), path("${sample_id}_dragonflye_${assembly_mode}_provenance.yml"),                             emit: provenance

    script:
      short_reads       = assembly_mode == "hybrid" ? "--R1 ${reads[0]} --R2 ${reads[1]}" : ""
      hybrid_long_reads = assembly_mode == "hybrid" ? "--reads ${reads[2]}" : ""
      long_reads        = assembly_mode == "long"   ? "--reads ${reads[0]}" : ""
      polypolish        = assembly_mode == "hybrid" ? "--polypolish" : ""
      """
      printf -- "- process_name: dragonflye\\n"                                                 >> ${sample_id}_dragonflye_${assembly_mode}_provenance.yml
      printf -- "  tools:\\n"                                                                   >> ${sample_id}_dragonflye_${assembly_mode}_provenance.yml
      printf -- "    - tool_name: dragonflye\\n"                                                >> ${sample_id}_dragonflye_${assembly_mode}_provenance.yml
      printf -- "      tool_version: \$(dragonflye --version | cut -d ' ' -f 2 | tr -d 'v')\\n" >> ${sample_id}_dragonflye_${assembly_mode}_provenance.yml

      dragonflye --cpus ${task.cpus} \
        --tmpdir './tmp' \
        ${short_reads} \
        ${hybrid_long_reads} \
	${long_reads} \
        --outdir ${sample_id}_assembly

      sed 's/^>/>${sample_id}_/' ${sample_id}_assembly/contigs.fa > ${sample_id}_dragonflye_${assembly_mode}.fa
      cp ${sample_id}_assembly/flye-unpolished.gfa ${sample_id}_dragonflye_${assembly_mode}_unpolished.gfa
      cp ${sample_id}_assembly/dragonflye.log ${sample_id}_dragonflye_${assembly_mode}.log
      """
}
