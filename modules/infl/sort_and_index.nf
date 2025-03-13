process sort_and_index {
    tag "sort_and_index:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), val(QC_status)

    output:
    tuple val(sampleId), path("${bam.baseName}_sorted.bam"), path("${bam.baseName}_sorted.bam.bai"), val(QC_status)

    script:
    def newBam = "${bam.baseName}_sorted.bam"
    """
    if [ ${QC_status} == "nie" ]; then
      touch "${bam.baseName}_sorted.bam"
      touch "${bam.baseName}_sorted.bam.bai"
    else
      samtools sort -o ${newBam} ${bam}
      samtools index ${newBam}
    fi
    """
}
