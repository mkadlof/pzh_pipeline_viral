process dehumanization_illumina  {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}*nohuman.fastq.gz"
    container  = params.main_image
    memory "20 GB"
    cpus 1
    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('dehumanized.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch nohuman.fastq.gz
      ERR_MSG="This module was eneterd with failed QC and poduced no valid output"
      parse_dehumanization.py --status "nie" --error "\${ERR_MSG}" -o dehumanized.json
    else
      samtools view mapped_reads.bam | cut -f1 | sort | uniq > lista_id_nohuman.txt
      # check if fastq file contains "\1" in R1, if so remove it, to match names of reads in bwa-produced bam file
       
      NUMBER_OF_READS=`zcat ${reads[0]} | grep '^@'`
      # Empty fastq will always have failed QC
      # however how many reads are checked must be considered dynamic 
      if [ \${NUMBER_OF_READS} -lt 100 ]; then
        TO_CHECK=`echo \${NUMBER_OF_READS} | awk '{print int(\$0 * 0.4)}'`
      else
        TO_CHECK=1000
      fi

      if [ \$(zcat ${reads[0]} | grep '^@' | head -n \${TO_CHECK} | awk 'BEGIN {i=0} { if (\$1 ~ /\\/1\$/ ) i+=1} END {print i}') -eq \${TO_CHECK} ] ; then
        seqkit replace -p "\\/1" -r "" ${reads[0]} | gzip > forward_paired_rename.fastq.gz
        seqkit replace -p "\\/2" -r "" ${reads[1]} | gzip > reverse_paired_rename.fastq.gz
  
        seqtk subseq forward_paired_rename.fastq.gz lista_id_nohuman.txt | gzip > ${sampleId}_forward_paired_nohuman.fastq.gz
        seqtk subseq reverse_paired_rename.fastq.gz lista_id_nohuman.txt | gzip > ${sampleId}_reverse_paired_nohuman.fastq.gz
        rm forward_paired_rename.fastq.gz reverse_paired_rename.fastq.gz
      else
        seqtk subseq ${reads[0]} lista_id_nohuman.txt | gzip > ${sampleId}_forward_paired_nohuman.fastq.gz
        seqtk subseq ${reads[1]} lista_id_nohuman.txt | gzip > ${sampleId}_reverse_paired_nohuman.fastq.gz
      fi

      find . -name "*paired_nohuman*gz" >> list_of_dehumanized_fastas.txt
      parse_dehumanization.py --status "tak" \
                              --input_fastas_list list_of_dehumanized_fastas.txt \
                              --output_path "${params.results_dir}/${sampleId}" \
                              --output dehumanized.json
    fi
    """
}

process dehumanization_nanopore {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}*nohuman.fastq.gz"
    container  = params.main_image
    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)
    cpus 1
    memory "20 GB"
    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('dehumanized.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch nohuman.fastq.gz
      ERR_MSG="This module was eneterd with failed QC and poduced no valid output"
      parse_dehumanization.py --status "nie" --error "\${ERR_MSG}" -o dehumanized.json
    else
      
      samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
      seqtk subseq ${reads} lista_id_nohuman.txt | gzip >> ${sampleId}_nohuman.fastq.gz
      find . -name "*nohuman.fastq.gz" >> list_of_dehumanized_fastas.txt
      parse_dehumanization.py --status "tak" \
                              --input_fastas_list list_of_dehumanized_fastas.txt \
                              --output_path "${params.results_dir}/${sampleId}" \
                              --output dehumanized.json

    fi
    """
}
