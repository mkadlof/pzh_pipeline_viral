process merging {
    tag "merging:${sampleId}"
    container  = params.main_image
    cpus { params.threads > 5 ? 5 : params.threads }
    memory "20 GB"
    input:
    tuple val(sampleId), path(filtering_bam), path(ivar_bam), val(QC_status)

    output:
    tuple val(sampleId), path('clean_sort_dedup_trimmed_sort.bam'), path('clean_sort_dedup_trimmed_sort.bam.bai'), val(QC_status)

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch clean_sort_dedup_trimmed_sort.bam
      touch clean_sort_dedup_trimmed_sort.bam.bai
    else
      samtools merge -o clean_sort_dedup_trimmed_sort_tmp.bam ${filtering_bam} ${ivar_bam}
      samtools sort -@ ${task.cpus} -o clean_sort_dedup_trimmed_sort.bam clean_sort_dedup_trimmed_sort_tmp.bam
      samtools index clean_sort_dedup_trimmed_sort.bam
    fi
    """
}
