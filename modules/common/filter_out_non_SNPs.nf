process filter_out_non_SNPs {
    // This simple modules keeps only SNPs in a vcf file
    tag "filtering out non-SNPs:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path('input.vcf.gz'), path('input.vcf.gz.tbi'), val(QC_status)

    output:
    tuple val(sampleId), path('output.vcf.gz'), path('output.vcf.gz.tbi'), emit: vcf

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch output.vcf.gz
      touch output.vcf.gz.tbi
    else
      zcat input.vcf.gz | awk '{if (\$1 ~ /^#/) {print \$0} else {if (length(\$5) == length(\$4)) {print \$0}}}' |  bcftools filter -O z -o output.vcf.gz -i "INFO/DP >= 1"

      tabix output.vcf.gz
    fi
    """
}
