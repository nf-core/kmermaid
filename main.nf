#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/kmermaid
========================================================================================
    Github : https://github.com/nf-core/kmermaid
    Website: https://nf-co.re/kmermaid
    Slack  : https://nfcore.slack.com/channels/kmermaid
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

include { KMERMAID } from './workflows/kmermaid'

//
// WORKFLOW: Run main nf-core/kmermaid analysis pipeline
//
workflow NFCORE_KMERMAID {
    KMERMAID ()
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
    NFCORE_KMERMAID ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
