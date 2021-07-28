#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/thatacseq
========================================================================================
    Github : https://github.com/nf-core/thatacseq
    Website: https://nf-co.re/thatacseq
    Slack  : https://nfcore.slack.com/channels/thatacseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { THATACSEQ } from './workflows/thatacseq'

//
// WORKFLOW: Run main nf-core/thatacseq analysis pipeline
//
workflow NFCORE_THATACSEQ {
    THATACSEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_THATACSEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
