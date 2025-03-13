process indelQual {
    tag "indelQual:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    // publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "forvariants.bam*"

    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(ref_genome)

    output:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai'), val(QC_status), emit: bam_and_qc
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai'), val(QC_status), path(ref_genome), emit: bam_genome_and_qc

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch forvariants.bam
      touch forvariants.bam.bai
    else
      lofreq indelqual --ref ${ref_genome} \
                       --out forvariants.bam \
                       --dindel ${bam}
      samtools sort -@ ${params.threads} \
                    -o forvariants.bam \
                    forvariants.bam
      samtools index forvariants.bam
    fi
    """
}
