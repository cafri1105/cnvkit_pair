#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.project_name = 'tcga+GRCh37_pair'
params.manifest      = '/mnt/NAS4/storage/ML_ECDNA/results/TCGA/WXS/manifest/final_TCGA_WXS+tumornormal_pair.csv'
params.bamDir        = '/mnt/NAS/storage/ML_ecDNA/tcga_wxs/GRCh37_bamfile'
params.norbamDir     = '/mnt/storage/ML_ecDNA/tcga_wxs_normal'
params.outputDir     = '/mnt/storage/ML_eCDNA/cnvkit_output'

// Entry Workflow
workflow {
    runCNVKitPair()
}

// Workflow Definition
workflow runCNVKitPair {

    main:
        println "Running CNVKit Pair Workflow!"

        // Step 1: Prepare Samples
        prepareSamples(params.manifest)

        // Step 2: Run CNVKit Analysis
        cnvkitAnalysis(prepareSamples.out)

    emit:
        out = cnvkitAnalysis.out
}

/*
================================================================================

                              Processes

================================================================================
*/

// -----------------------------------------------------------------
// prepareSamples
// -----------------------------------------------------------------
process prepareSamples {

    tag "Prepare samples from manifest"

    input:
        path manifest

    output:
        path "prepared_samples.tsv"

    script:
    """
    #!/bin/bash
    awk 'BEGIN {FS=","; OFS="\t"}
    NR == 1 {
        # Find the column indices for name_tumor, name_normal, and reference_genome
        for (i = 1; i <= NF; i++) {
            if (\$i == "name_tumor") tumor_idx = i
            if (\$i == "name_normal") normal_idx = i
            if (\$i == "reference_genome") genome_idx = i
        }
        # Always print the header row
        print "name_normal", "name_tumor", "reference_genome"
        next
    }
    {
        # Output the columns in the specified order
        print \$normal_idx, \$tumor_idx, \$genome_idx
    }' ${manifest} > prepared_samples.tsv
    """
}

// -----------------------------------------------------------------
// cnvkitAnalysis
// -----------------------------------------------------------------
process cnvkitAnalysis {

    tag "Run CNVKit Analysis"
    container "/mnt/NAS4/home/hj/sing_img/hj_cnvlit.sif"
    publishDir "${params.outputDir}/${params.project_name}/tumor_pair/cnvkit/", mode: 'copy'

    input:
        path "prepared_samples.tsv"

    output:
        path "cnvkit_results/"

    script:
    """
    #!/bin/bash
    echo "Running CNVKit analysis on prepared_samples.tsv"
    mkdir -p cnvkit_results

    # Read the input TSV file
    tail -n +2 prepared_samples.tsv | while IFS=\$'\t' read -r name_normal name_tumor reference_genome; do

        tumor_bamdir="${params.bamDir}"
        normal_bamdir="${params.norbamDir}"

        name_tumor_base=\$(basename "\${name_tumor}" .bam)
        name_normal_base=\$(basename "\${name_normal}" .bam)
        outdir="cnvkit_results/\${name_tumor_base}_\${name_normal_base}"

        if [[ -f "/mnt/storage/ML_ecDNA/cns_file/\${name_tumor_base}_\${name_normal_base}.call.cns" ]] || \
           [[ -f "\${outdir}/\${name_tumor_base}_\${name_normal_base}.call.cns" ]]; then
            echo "Skipping \${name_tumor} as .call.cns file already exists."
            continue
        fi

        if [[ "\${reference_genome}" =~ ^(HG18_Broad_variant|NCBI36_BCM_variant|NCBI36_WUGSC_variant)\$ ]]; then
            echo "Skipping \${name_tumor} due to unsupported genome version: \${reference_genome}"
            continue
        fi

        case "\${reference_genome}" in
            HG19_Broad_variant)
                GENOME="/mnt/NAS/storage/references/GRCh37/HG19_Broad_variant.fasta"
                ;;
            GRCh37-lite_WUGSC_variant_1)
                GENOME="/mnt/NAS/storage/references/GRCh37/GRCh37-lite_WUGSC_variant_1.fa.gz"
                ;;
            GRCh37-lite)
                GENOME="/mnt/NAS/storage/references/GRCh37/GRCh37-lite.fa"
                ;;
            GRCh37-lite-+-HPV_Redux-build)
                GENOME="/mnt/NAS/storage/references/GRCh37/GRCh37-lite-+-HPV_Redux-build.fa.gz"
                ;;
            *)
                echo "Error: Unknown reference genome: \${reference_genome}"
                continue
                ;;
        esac

        mkdir -p "\${outdir}"
        genome_bed_path="\${outdir}/genome.bed"
        target_bed="\${outdir}/targets.bed"
        antitarget_bed="\${outdir}/antitargets.bed"

        cnvkit.py access "\${GENOME}" -o "\${genome_bed_path}"
        cnvkit.py target /mnt/NAS4/home/hj/reference/cnvkit_data/cleaned_output_file.interval_list --split -o "\${target_bed}"
        cnvkit.py antitarget /mnt/NAS4/home/hj/reference/cnvkit_data/cleaned_output_file.interval_list -g "\${genome_bed_path}" -o "\${antitarget_bed}"

        cnvkit.py coverage "\${tumor_bamdir}/\${name_tumor}" "\${target_bed}" -o "\${outdir}/\${name_tumor_base}.targetcoverage.cnn"
        cnvkit.py coverage "\${tumor_bamdir}/\${name_tumor}" "\${antitarget_bed}" -o "\${outdir}/\${name_tumor_base}.antitargetcoverage.cnn"
        cnvkit.py coverage "\${normal_bamdir}/\${name_normal}" "\${target_bed}" -o "\${outdir}/\${name_normal_base}.targetcoverage.cnn"
        cnvkit.py coverage "\${normal_bamdir}/\${name_normal}" "\${antitarget_bed}" -o "\${outdir}/\${name_normal_base}.antitargetcoverage.cnn"

        cnvkit.py batch -p 8 "\${tumor_bamdir}/\${name_tumor}" \
            --normal "\${normal_bamdir}/\${name_normal}" \
            --targets "\${target_bed}" \
            --fasta "\${GENOME}" \
            --access "\${genome_bed_path}" \
            --drop-low-coverage \
            --output-reference "\${outdir}/\${name_tumor_base}_\${name_normal_base}.reference.cnn"

        cnvkit.py fix "\${outdir}/\${name_tumor_base}.targetcoverage.cnn" \
            "\${outdir}/\${name_tumor_base}.antitargetcoverage.cnn" \
            "\${outdir}/\${name_tumor_base}_\${name_normal_base}.reference.cnn" \
            -o "\${outdir}/\${name_tumor_base}_\${name_normal_base}.cnr"

        cnvkit.py segment "\${outdir}/\${name_tumor_base}_\${name_normal_base}.cnr" \
            --method cbs \
            --drop-low-coverage \
            --min-variant-depth 20 \
            -o "\${outdir}/\${name_tumor_base}_\${name_normal_base}.cns"

        cnvkit.py call "\${outdir}/\${name_tumor_base}_\${name_normal_base}.cns" \
            -o "\${outdir}/\${name_tumor_base}_\${name_normal_base}.call.cns"

        cnvkit.py export seg "\${outdir}/\${name_tumor_base}_\${name_normal_base}.call.cns" \
            -o "\${outdir}/\${name_tumor_base}_\${name_normal_base}.seg"

    done
    """
}
