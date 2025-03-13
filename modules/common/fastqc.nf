process fastqc {
    tag "${prefix}_fastqc:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/QC", mode: 'copy', pattern: "*_reads_quality_histogram.csv"
    publishDir "${params.results_dir}/${sampleId}/QC", mode: 'copy', pattern: "*_reads_length_histogram.csv"
    publishDir "${params.results_dir}/${sampleId}/QC", mode: 'copy', pattern: "*_position_quality_plot.csv"
    cpus { params.threads > 10 ? 10 : params.threads }
    memory '20 GB'
    // publishDir "${params.results_dir}/${sampleId}/json_output", mode: 'copy', pattern: "*.json"
    container  = params.main_image

    input:
    tuple val(sampleId), path(reads), val(QC_STATUS)
    val(prefix)

    output:
    tuple val(sampleId), path("*csv"),                                                   emit: publishdir // Wykresy, kopiujemy bo json wskazuje do publishdir
    tuple val(sampleId), path("forward_${prefix}.json"), path("reverse_${prefix}.json"), emit: json // Sam json do kopiowania w celu zlozenia "ostatecznego jsona"
    tuple val(sampleId), env(QC_STATUS_EXIT),                                            emit: qcstatus // Sam QC status

    script:
    if (QC_STATUS == null) { QC_STATUS="tak" } // Domyslna wartosc w przypadku gdy user nie poda QC status
    """
    # Set up QC_STATUS to "tak" if a given module does not provide this value via input
    # run_fastqc_and_generate_json.py will always produce output required by these module, even if no valid fastq file is provided, or fastq file
    # does not meet predeifned criteria. The script returns to values status (tak, nie, blad) and total numberof bases in fastqfile (0 if status is not tak)
    if [ ${QC_STATUS} == "nie" ]; then
      ERROR_MSG="Initial QC received by this module was nie"
      touch dummy.csv
    else
      ERROR_MSG=""
    fi

    DANE_FORWARD=(`run_fastqc_and_generate_json.py -i ${reads[0]} -m 4048 -c ${task.cpus} -x ${params.min_number_of_reads} -y ${params.min_median_quality} -s ${QC_STATUS} -r "\${ERROR_MSG}" -e ${prefix} -p "${params.results_dir}/${sampleId}/QC" -o forward_${prefix}.json`)
    STATUS_FORWARD_ALL="\${DANE_FORWARD[0]}"
    BASES_FORWARD="\${DANE_FORWARD[1]}"
    DANE_REVERSE=(`run_fastqc_and_generate_json.py -i ${reads[1]} -m 4048 -c ${task.cpus} -x ${params.min_number_of_reads} -y ${params.min_median_quality} -s ${QC_STATUS} -r "\${ERROR_MSG}" -e ${prefix} -p "${params.results_dir}/${sampleId}/QC" -o reverse_${prefix}.json`)
    STATUS_REVERSE_ALL="\${DANE_REVERSE[0]}"
    BASES_REVERSE="\${DANE_REVERSE[1]}"
    TOTAL_BASES=`echo "\${BASES_FORWARD} + \${BASES_REVERSE}" | bc -l`
    if [[ \${STATUS_FORWARD_ALL} == "nie"  || \${STATUS_REVERSE_ALL} == "nie"  || \${STATUS_FORWARD_ALL} == "blad"  || \${STATUS_REVERSE_ALL} == "blad" ]]; then
      QC_STATUS_EXIT="nie" # moduly "nizej" dostaja status nie
    else
      QC_STATUS_EXIT="tak"
    fi
    """
}


process run_fastqc_nanopore {
  // Modification of run_fastqc_illumina
  // That requires one fastq.gz filei
  
  tag "fastqc for sample ${sampleId}"
  container  = params.main_image
  cpus { params.threads > 10 ? 10 : params.threads }
  memory '20 GB'
  publishDir "${params.results_dir}/${sampleId}/QC", mode: 'copy', pattern: "*reads_quality_histogram.csv"
  publishDir "${params.results_dir}/${sampleId}/QC", mode: 'copy', pattern: "*reads_length_histogram.csv"
  publishDir "${params.results_dir}/${sampleId}/QC", mode: 'copy', pattern: "*position_quality_plot.csv"
  // publishDir "${params.results_dir}/${sampleId}/json_output", mode: 'copy', pattern: "*.json"

  input:
  tuple val(sampleId), path(reads), val(QC_STATUS)
  val(prefix)

  output:
  tuple val(sampleId), path("*csv"), emit: publishdir
  tuple val(sampleId), path("forward.json"), emit: json
  tuple val(sampleId), env(QC_STATUS_EXIT), emit: qcstatus
  tuple val(sampleId), env(QC_STATUS_EXIT), env(TOTAL_BASES), emit: qcstatus_and_values

  script:
  if (QC_STATUS == null) { QC_STATUS="tak" } //
  """
  if [ ${QC_STATUS} == "nie" ]; then
    ERROR_MSG="Initial QC received by this module was nie"
  else
    ERROR_MSG=""
  fi
 
  DANE_FORWARD=(`run_fastqc_and_generate_json.py -i ${reads} -m 4048 -c ${task.cpus} -x ${params.min_number_of_reads} -y ${params.min_median_quality} -s ${QC_STATUS} -r "\${ERROR_MSG}" -e ${prefix} -p "${params.results_dir}/${sampleId}/QC" -o forward.json`)
  STATUS_FORWARD="\${DANE_FORWARD[0]}"
  TOTAL_BASES="\${DANE_FORWARD[1]}"

  if [ \${STATUS_FORWARD} != "tak" ]; then
    QC_STATUS_EXIT="nie"
    touch reads_quality_histogram.csv
    touch reads_length_histogram.csv
    touch position_quality_plot.csv
  else
    QC_STATUS_EXIT="tak"
  fi

  """
}
