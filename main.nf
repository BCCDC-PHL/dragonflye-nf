#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files }                 from './modules/hash_files.nf'
include { fastp }                      from './modules/fastp.nf'
include { fastp_json_to_csv }          from './modules/fastp.nf'
include { filtlong }                   from './modules/long_read_qc.nf'
include { nanoq as nanoq_pre_filter }  from './modules/long_read_qc.nf'
include { nanoq as nanoq_post_filter } from './modules/long_read_qc.nf'
include { merge_nanoq_reports }        from './modules/long_read_qc.nf'
include { dragonflye }                 from './modules/dragonflye.nf'
include { prokka }                     from './modules/prokka.nf'
include { bakta }                      from './modules/bakta.nf'
include { quast }                      from './modules/quast.nf'
include { parse_quast_report }         from './modules/quast.nf'
include { bandage }                    from './modules/long_read_qc.nf'
include { pipeline_provenance }        from './modules/provenance.nf'
include { collect_provenance }         from './modules/provenance.nf'


workflow {
  ch_start_time = Channel.of(LocalDateTime.now())
  ch_pipeline_name = Channel.of(workflow.manifest.name)
  ch_pipeline_version = Channel.of(workflow.manifest.version)

  ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

  if (params.samplesheet_input != 'NO_FILE') {
    if (params.hybrid) {
      ch_assembly_mode = Channel.of("hybrid")
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2'], it['LONG']]] }
      ch_short_reads = ch_fastq.map{ it -> [it[0], [it[1][0], it[1][1]]] }
      ch_long_reads = ch_fastq.map{ it -> [it[0], it[1][2]] }
    } else if (params.long_only) {
      ch_assembly_mode = Channel.of("long")
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['LONG']]] }
      ch_short_reads = Channel.of()
      ch_long_reads = ch_fastq
    } else {
      ch_assembly_mode = Channel.of("short")
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2']]] }
      ch_short_reads = ch_fastq
      ch_long_reads = Channel.of()
    }
  } else {
    if (params.hybrid) {
      ch_assembly_mode = Channel.of("hybrid")
      ch_short_reads = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }.unique{ it -> it[0] }
      ch_long_reads = Channel.fromPath( params.long_reads_search_path ).map{ it -> [it.baseName.split("_")[0], [it]] }
      ch_fastq = ch_short_reads.join(ch_long_reads).map{ it -> [it[0], it[1] + it[2]] }
    } else if (params.long_only) {
      ch_assembly_mode = Channel.of("long")
      ch_short_reads = Channel.of()
      ch_long_reads = Channel.fromPath( params.long_reads_search_path ).map{ it -> [it.baseName.split("_")[0], [it]] }
      ch_fastq = ch_long_reads
    } else {
      ch_assembly_mode = Channel.of("short")
      ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }.unique{ it -> it[0] }
      ch_short_reads = ch_fastq
      ch_long_reads = Channel.of()
    }
  }

  if (params.long_only && params.hybrid) {
    System.out.println("Choose one of --long or --hybrid (but not both).")
    System.exit(-1)
  }  

  main:
    ch_provenance = ch_fastq.map{ it -> it[0] }

    hash_files(ch_fastq.combine(Channel.of("fastq-input")))

    if (params.hybrid) {
      fastp(ch_short_reads)
      fastp_json_to_csv(fastp.out.json)
      nanoq_pre_filter(ch_long_reads.combine(Channel.of("pre_filter")))
      filtlong(ch_long_reads)
      nanoq_post_filter(filtlong.out.filtered_reads.combine(Channel.of("post_filter")))
      merge_nanoq_reports(nanoq_pre_filter.out.report.join(nanoq_post_filter.out.report))
      dragonflye(fastp.out.trimmed_reads.join(filtlong.out.filtered_reads).map{ it -> [it[0], [it[1], it[2], it[3]]] }.combine(ch_assembly_mode))
    } else if (params.long_only) {
      nanoq_pre_filter(ch_long_reads.combine(Channel.of("pre_filter")).map{ it -> [it[0], it[1][0], it[2]] })
      filtlong(ch_long_reads)
      nanoq_post_filter(filtlong.out.filtered_reads.combine(Channel.of("post_filter")))
      merge_nanoq_reports(nanoq_pre_filter.out.report.join(nanoq_post_filter.out.report))
      dragonflye(ch_long_reads.combine(ch_assembly_mode))
    }

    if (params.prokka) {
      prokka(dragonflye.out.assembly)
    }

    if (params.bakta) {
      bakta(dragonflye.out.assembly)
    }

    quast(dragonflye.out.assembly)
    bandage(dragonflye.out.assembly_graph)

    parse_quast_report(quast.out.tsv)

    //
    // Provenance collection processes
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    if (params.hybrid || params.long_only) {
      ch_provenance = ch_provenance.join(nanoq_pre_filter.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
      ch_provenance = ch_provenance.join(filtlong.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
      ch_provenance = ch_provenance.join(nanoq_post_filter.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    if (!params.long_only) {
      ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    ch_provenance = ch_provenance.join(dragonflye.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    if (params.prokka) {
      ch_provenance = ch_provenance.join(prokka.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    if (params.bakta) {
      ch_provenance = ch_provenance.join(bakta.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    ch_provenance = ch_provenance.join(quast.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
