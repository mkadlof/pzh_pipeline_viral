process picard_wgsMetrics {
    tag "wgsMetrics:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "${sampleId}_coverage_barplot_*"
    publishDir "${params.results_dir}/${sampleId}/", mode: 'copy', pattern: "${sampleId}_coverage_histogram.csv"

    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(ref_genome), path('Primer_usage.txt')

    output:
    tuple val(sampleId), path("viral_genome_data.json"), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      ERR_MSG="QC failed: an error occurred in a prior processing step"
      picard_parser.py --status "nie" \
                       --output viral_genome_data.json \
                       --error "\${ERR_MSG}"
    else
      java -jar /opt/picard/picard.jar CollectWgsMetrics --REFERENCE_SEQUENCE ${ref_genome} \
                                                         --MINIMUM_BASE_QUALITY ${params.quality_initial} \
                                                         --MINIMUM_MAPPING_QUALITY ${params.min_mapq} \
                                                         --INPUT ${bam} \
                                                         --OUTPUT picard_statistics.txt

      bedtools genomecov -d -ibam ${bam} > genomecov.bedgraph

      # Split the genomecov.bedgraph file into segments using awk
     
      input_file="genomecov.bedgraph"
      output_prefix="${sampleId}_coverage_barplot_"
      summary_coverage_file="coverage_barplot_files.txt"

      awk 'BEGIN { OFS = "," }
      {
      oryg=\$0;
      split(oryg, a, "\\t")
      gsub("/", "_");
      gsub("\\\\.", "_");
      if (\$1 != prev) {
        if (prev != "") {
          close(out)
         }
         prev = \$1
         out = "'"\${output_prefix}"'" prev ".csv"
         out_summary = "'"\${summary_coverage_file}"'"
         i = 1
         print "#indeks", "segment", "pozycja", "pokrycie" >> out
         print a[0], out >> out_summary
         
      }
      print i, a[1], a[2], a[3] >> out;
      i += 1
      }' "\${input_file}"


      rm genomecov.bedgraph

      picard_parser.py --input_file_picard picard_statistics.txt \
                      --input_file_primers Primer_usage.txt \
                      --input_file_bedgraph \${summary_coverage_file} \
                      --output_path "${params.results_dir}/${sampleId}/" \
                      --status "tak" \
                      --output viral_genome_data.json


    fi
    """
}
