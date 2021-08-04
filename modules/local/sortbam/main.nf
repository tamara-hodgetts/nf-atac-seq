// Import generic module functions
// include { initOptions; saveFiles; getSoftwareName } from './functions'

// params.options = [:]
// options        = initOptions(params.options)

process SORT_BAM {
        tag "$name"
        label 'process_medium'
        if (params.save_align_intermeds) {
            publishDir path: "${params.outdir}/bwa/library", mode: params.publish_dir_mode,
                saveAs: { filename ->
                          if (filename.endsWith('.flagstat')) "samtools_stats/$filename"
                          else if (filename.endsWith('.idxstats')) "samtools_stats/$filename"
                          else if (filename.endsWith('.stats')) "samtools_stats/$filename"
                          else filename
                        }
        }
    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0"
    }
    
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path('*.sorted.{bam,bam.bai}') into ch_sort_bam_merge
    path '*.{flagstat,idxstats,stats}' into ch_sort_bam_flagstat_mqc

    script:
    prefix = "${name}.Lb"
    """
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    """
    }