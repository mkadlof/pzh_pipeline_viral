process coinfection_analysis {
    tag "coinfection_analysis:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "*allele_usage_histogram.txt"

    input:
    tuple val(sampleId), path(detected_variants_varscan_coinfection_txt), val(QC_status)

    output:
    tuple val(sampleId), path("${sampleId}_allele_usage_histogram.txt"), emit: to_pubdir
    tuple val(sampleId), path('custom_coinfection_analysis.json'), emit: json

    script:
    """
    # The coinfection analysis is based on similarity to 3 known samples from EQA2023.
    # The new version of the script allows specifying any number of samples that are known to be
    # coinfected. Simply add additional files as positional arguments.

    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_allele_usage_histogram.txt
      coinfection_status="nie"
      ERR_MSG="This module was eneterd with failed QC and poduced no valid output"
      echo -e "{\\"coinfection_status\\":\\"\${coinfection_status}\\", \
                \\"coinfetion_error_message\\":\\"\${ERR_MSG}\\"}" >> custom_coinfection_analysis.json
    else
      # use different files for nanopore and illumina data
      if [ ${params.machine} == 'Illumina' ]; then
        RESULTS=(`predict_coinfection_illumina.py ${detected_variants_varscan_coinfection_txt} \
                                        ${sampleId} \
                                        /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.09_coinfections.txt \
                                        /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.17_coinfections.txt \
                                        /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS2.32_coinfections.txt`)
      elif [ ${params.machine} == 'Nanopore' ]; then
        RESULTS=(`predict_coinfection_illumina.py ${detected_variants_varscan_coinfection_txt} \
                                        ${sampleId} \
                                        /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS1.04_coinfections.txt \
                                        /home/data/sarscov2/coinfections/ESIB_EQA_2023.SARS1.05_coinfections.txt`)
      fi

      cp allele_usage_histogram.txt ${sampleId}_allele_usage_histogram.txt
      coinfection_status="tak"
      coinfection_result=\${RESULTS[0]}
      coinfection_pvalue=`echo \${RESULTS[1]} | awk '{printf "%.2f", \$0}'`
      coinfection_histogram_file="${params.results_dir}/${sampleId}/${sampleId}_allele_usage_histogram.txt"

      echo -e "{\\"coinfection_status\\":\\"\${coinfection_status}\\",
                \\"coinfection_result\\":\\"\${coinfection_result}\\",
                \\"coinfection_pvalue\\":\${coinfection_pvalue},
                \\"coinfection_histogram_file\\":\\"\${coinfection_histogram_file}\\"}" >> custom_coinfection_analysis.json

    fi

    """
}
