process bwa {
    tag "bwa:${sampleId}"
    cpus params.threads
    container  = params.main_image
    memory "20 GB"
    input:
    tuple val(sampleId), path(reads), path(ref_genome_with_index), path("primers.bed"), val(QC_status)

    output:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), env(QC_exit), emit: only_bam
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path(ref_genome_with_index), env(QC_exit), emit: bam_and_genome
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), path(ref_genome_with_index), path("primers.bed"), env(QC_exit), emit: to_coinfection
    tuple val(sampleId), path("mapping.json"), emit: json
    script:
    // Check the index of a file with fasta extension in ref_genome_with_index list
    
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }

    """
    if [ ${QC_status} == "nie" ]; then
      touch mapped_reads.bam
      touch mapped_reads.bam.bai
      QC_exit="nie"
      ERR_MSG="This module was eneterd with failed QC and poduced no valid output"

      echo -e "{\\"status\\":\\"\${QC_exit}\\", \
                \\"error_message\\": \\"\${ERR_MSG}\\"}" >> mapping.json
    else
      if [ final_index -eq -1 ]; then
        touch mapped_reads.bam
        touch mapped_reads.bam.bai
        QC_exit="nie"
        ERR_MSG="Upstream module did not provide a valid genome file"
        echo -e "{\\"status\\":\\"\${QC_exit}\\", \
                  \\"error_message\\": \\"\${ERR_MSG}\\"}" >> mapping.json
      else
        bwa mem -t ${params.threads} -T ${params.min_mapq} ${ref_genome_with_index[final_index]} ${reads[0]} ${reads[1]} | \
        samtools view -@ ${params.threads} -Sb -f 3 -F 2048 - | \
        samtools sort -@ ${params.threads} -o mapped_reads.bam -
        samtools index mapped_reads.bam
        NO_READS=`samtools view mapped_reads.bam | wc -l`
        if [ \${NO_READS} -lt ${params.min_number_of_reads} ]; then
          QC_exit="nie"
          ERR_MSG="There are less than ${params.min_number_of_reads} mapping to reference genome"
          echo -e "{\\"status\\":\\"blad\\", \
                    \\"error_message\\": \\"\${ERR_MSG}\\"}" >> mapping.json
        else
          QC_exit="tak"
          echo -e "{\\"status\\":\\"\${QC_exit}\\", \
                    \\"mapped_reads_number\\": \${NO_READS}}" >> mapping.json
        fi  # if na brak poprawnych odczytow po mapowaniu
      fi # if na nie podanie poprawnego genomu
    fi # if na QC status
    """
}
