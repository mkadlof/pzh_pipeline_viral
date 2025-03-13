process minimap2 {
    container  = params.main_image 
    tag "minimap2:${sampleId}"
    cpus params.threads
    memory "20 GB"
    input:
    tuple val(sampleId), path(reads), path("ref_genome.fasta"), path("primers.bed"), val(QC_status)

    output:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), env(QC_exit), emit: bam_and_qc
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path("ref_genome.fasta"), env(QC_exit), emit: bam_and_genome
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), env(QC_exit), path("ref_genome.fasta"), emit: for_snpeff
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path("ref_genome.fasta"), path("primers.bed"), env(QC_exit), emit: bam_and_genome_and_primers
    tuple val(sampleId), path("mapping.json"), emit: json

    script:

    """
    if [ ${QC_status} == "nie" ]; then
      touch mapped_reads.bam
      touch mapped_reads.bam.bai
      QC_exit="nie"
      ERR_MSG="This module was eneterd with failed QC and poduced no valid output"
      echo -e "{\\"status\\":\\"\${QC_exit}\\", \
                \\"error_message\\": \\"\${ERR_MSG}\\"}" >> mapping.json
     
   else
      minimap2 -a -x map-ont -t ${params.threads} -o tmp.sam ref_genome.fasta ${reads}
      samtools view -@ ${params.threads} -Sb -F 2052 tmp.sam | \
      samtools sort -@ ${params.threads} -o mapped_reads.bam -
      samtools index mapped_reads.bam
      rm tmp.sam
      NO_READS=`samtools view mapped_reads.bam | wc -l`
      if [ \${NO_READS} -lt ${params.min_number_of_reads} ]; then
        QC_exit="nie"
        ERR_MSG="There are less than ${params.min_number_of_reads} mapping to reference genome"
        echo -e "{\\"status\\":\\"blad\\", \
                  \\"error_message\\": \\"\${ERR_MSG}\\"}" >> mapping.json
      else
        QC_exit="tak"
        echo -e "{\\"status\\":\\"tak\\", \
                  \\"mapped_reads_number\\": \${NO_READS}}" >> mapping.json
      fi  # if na brak poprawnych odczytow po mapowaniu
    fi # if na QC status
    """
}
