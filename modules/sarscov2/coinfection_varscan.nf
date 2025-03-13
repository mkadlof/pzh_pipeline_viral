process coinfection_varscan {
    tag "coinfection_varscan:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path(for_contamination_mpileup), val(QC_status)

    output:
    tuple val(sampleId), path('detected_variants_varscan_contamination.txt'), val(QC_status), emit: to_pubdir
    

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch detected_variants_varscan_contamination.txt
    else
      if [ ${params.machine} == 'Illumina' ]; then
        min_cov=20
        varscan_qual=30
      elif [ ${params.machine} == 'Nanopore' ]; then
        min_cov=50
        varscan_qual=1
      fi
      java -jar /opt/varscan/VarScan.v2.4.6.jar pileup2snp ${for_contamination_mpileup} \
                                                     --min-avg-qual \${varscan_qual} \
                                                     --p-value 0.9 \
                                                     --min-var-freq 0.05 \
                                                     --min-coverage \${min_cov} \
                                                     --variants \
                                                     --min-reads2 0 > detected_variants_varscan_contamination.txt
    fi
    """
}
