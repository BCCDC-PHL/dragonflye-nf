process rotate_replicons {

    tag { sample_id + ' / ' + assembly_mode }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_${assembler}_${assembly_mode}_rotated.fa", mode: 'copy'
    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_${assembler}_${assembly_mode}_replicon_rotation.csv", mode: 'copy'

    input:
      tuple val(sample_id), path(assembly), val(assembler), val(assembly_mode)

    output:
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_rotated.fa"),  val("dragonflye"), val(assembly_mode), emit: rotated_assembly
      tuple val(sample_id), path("${sample_id}_${assembler}_${assembly_mode}_replicon_rotation.csv"), emit: replicon_rotation

    script:
    """
    mkdir -p ./tmp

    rotate-replicons.py \
      --blast-threads ${task.cpus} \
      --input ${assembly} \
      --start-genes ${params.start_genes} \
      --tmp ./tmp \
      --no-cleanup \
      --output ${sample_id}_${assembler}_${assembly_mode}_rotated.fa \
      > ${sample_id}_${assembler}_${assembly_mode}_replicon_rotation.csv
    """
}