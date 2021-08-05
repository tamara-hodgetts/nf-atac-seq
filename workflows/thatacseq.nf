/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowThatacseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

//ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
 
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

// def multiqc_options   = modules['multiqc']
// multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
// include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

// including additional modules
include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'
include { BWA_INDEX } from '../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM } from '../modules/nf-core/modules/bwa/mem/main'

include { SAMTOOLS_SORT } from '../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_FLAGSTAT } from '../modules/nf-core/modules/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../modules/nf-core/modules/samtools/idxstats/main'
include { SAMTOOLS_STATS } from '../modules/nf-core/modules/samtools/stats/main'

//include {samtools_index; samtools_view; samtools_faidx; samtools_sort} from '../modules/nf-core/modules/samtools'


// building a BWA index
// if (!params.bwa_index) {
//     process BWA_INDEX {
//         tag "$fasta"
//         label 'process_high'
//         publishDir path: { params.save_reference ? "${params.outdir}/genome" : params.outdir },
//             saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

//         input:
//         path fasta

//         output:
//         path 'BWAIndex' into ch_bwa_index

//         script:
//         """
//         bwa index -a bwtsw $fasta
//         mkdir BWAIndex && mv ${fasta}* BWAIndex
//         """
//     }
// }
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow THATACSEQ {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split("_")[0..-2].join("_")
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .set { ch_fastq }
    //ch_fastq | view

    //
    // MODULE: Run FastQC
    //
     FASTQC (
         INPUT_CHECK.out.reads
     )
      ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    //ch_software_versions | view
    //
    TRIMGALORE (
        INPUT_CHECK.out.reads
    )
    //
    BWA_INDEX (
        params.fasta
    )
    //
     BWA_MEM (
          INPUT_CHECK.out.reads, BWA_INDEX.out.index
    )
    // BWA_MEM.out.bam | view
    //
    //
    SAMTOOLS_SORT (
         BWA_MEM.out.bam
    )
    //
    SAMTOOLS_INDEX (
        BWA_MEM.out.bam
    )
    // 
    SAMTOOLS_FLAGSTAT (
        BWA_MEM.out.bam, SAMTOOLS_INDEX.out.bai
    )
    //
    SAMTOOLS_IDXSTATS (
        BWA_MEM.out.bam, SAMTOOLS_INDEX.out.bai
    )
    // 
    SAMTOOLS_STATS (
        BWA_MEM.out.bam, SAMTOOLS_INDEX.out.bai
    )
    // converting bam file to a sorted BAM
    // process SORT_BAM {
    //     tag "$name"
    //     label 'process_medium'
    //     if (params.save_align_intermeds) {
    //         publishDir path: "${params.outdir}/bwa/library", mode: params.publish_dir_mode,
    //             saveAs: { filename ->
    //                       if (filename.endsWith('.flagstat')) "samtools_stats/$filename"
    //                       else if (filename.endsWith('.idxstats')) "samtools_stats/$filename"
    //                       else if (filename.endsWith('.stats')) "samtools_stats/$filename"
    //                       else filename
    //                     }
    //     }

    // input:
    // tuple val(name), path(bam)

    // output:
    // tuple val(name), path('*.sorted.{bam,bam.bai}') into ch_sort_bam_merge
    // path '*.{flagstat,idxstats,stats}' into ch_sort_bam_flagstat_mqc

    // script:
    // prefix = "${name}.Lb"
    // """
    // samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
    // samtools index ${prefix}.sorted.bam
    // samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    // samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    // samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    // """
    // }

    // 
    // //
    // // MODULE: Pipeline reporting
    // //
    // ch_software_versions
    //     .map { it -> if (it) [ it.baseName, it ] }
    //     .groupTuple()
    //     .map { it[1][0] }
    //     .flatten()
    //     .collect()
    //     .set { ch_software_versions }

    // GET_SOFTWARE_VERSIONS (
    //     ch_software_versions.map { it }.collect()
    // )

    // //
    // // MODULE: MultiQC
    // //
    // workflow_summary    = WorkflowThatacseq.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report       = MULTIQC.out.report.toList()
    // ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
