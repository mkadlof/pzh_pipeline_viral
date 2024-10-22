process run_fastqc_nanopore {
  // Modification of run_fastqc_illumina
  // That requires one fastq.gz file
  tag "fastqc for sample ${x}"
  publishDir "${params.results_dir}/${sampleId}/QC", mode: 'copy'
  publishDir "${params.results_dir}/${sampleId}/json_output", mode: 'copy', pattern: "*.json"

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
 
  DANE_FORWARD=(`run_fastqc_and_generate_json.py -i ${reads} -m ${params.memory} -c ${params.threads} -x ${params.min_number_of_reads} -y ${params.min_median_quality} -s ${QC_STATUS} -r "\${ERROR_MSG}" -e ${prefix} -p "${params.results_dir}/${sampleId}/QC" -o forward.json`)
  STATUS_FORWARD="\${DANE_FORWARD[0]}"
  TOTAL_BASES="\${DANE_FORWARD[1]}"

  if [ \${STATUS_FORWARD} != "tak" ]; then
    QC_STATUS_EXIT="nie"
  else
    QC_STATUS_EXIT="tak"
  fi

  """
}