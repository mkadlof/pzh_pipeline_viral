process medaka_first_round {
    // Odplaenie medaki podczas pierwszej rundy analizy
    tag "medaka_first_round:${sampleId}"
    container  = params.medaka_image
    cpus params.threads
    memory "20 GB"
    input:
    tuple val(sampleId), path('trimmed.bam'), path('trimmed.bam.bai'), val(QC_status), path('genome.fasta')
    output:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz'), path('medaka_annotated_filtered.vcf.gz.tbi'),  path('genome.fasta'), val(QC_status), emit: vcf

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch medaka_annotated_filtered.vcf.gz
      touch medaka_annotated_filtered.vcf.gz.tbi
    else
      medaka inference --model ${params.medaka_model} \
                       --threads ${params.threads} \
                       --chunk_len ${params.medaka_chunk_len} \
                       --chunk_ovlp ${params.medaka_chunk_overlap} \
                       trimmed.bam forvariants.hdf
     
      medaka vcf forvariants.hdf genome.fasta medaka_initial.vcf
      medaka tools annotate medaka_initial.vcf genome.fasta trimmed.bam medaka_annotated.vcf

      bgzip medaka_annotated.vcf; tabix medaka_annotated.vcf.gz
      QUAL=`echo ${params.first_round_pval} | awk '{print 10*-log(\$1)/log(10)}'` # dodanie int powodowalo blad w zaookraglaniu, teraz dla 0.1 program poprawnie zwoci wartosc 10

      bcftools filter -O z -o medaka_annotated_filtered.vcf.gz -i "GQ >= \${QUAL} && DP >= ${params.min_cov}" medaka_annotated.vcf.gz
      tabix medaka_annotated_filtered.vcf.gz
    fi
    """
}

process medaka_second_round {
    // W drugiej rundzie jestesmy dokladniejsi i stosujemy analize ODDZIELNIE
    // Dla kazdego segmentu genomu
    // Co pozwala dobrac inne wartosci --chunk_len i --chunk_ovlp
    tag "medaka_second_roud:${sampleId}"
    container  = params.medaka_image
    cpus params.threads
    memory "20 GB"
    input:
    tuple val(sampleId), path('trimmed.bam'), path('trimmed.bam.bai'), val(QC_status), path('genome.fasta'), path('pre_trimmed.bam'), path('pre_trimmed.bam.bai')
    output:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz'), path('medaka_annotated_filtered.vcf.gz.tbi'), path('medaka_annotated.vcf.gz'), path('medaka_annotated.vcf.gz.tbi'), path('genome.fasta'), val(QC_status), emit: vcf

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch medaka_annotated_filtered.vcf.gz
      touch medaka_annotated_filtered.vcf.gz.tbi
      touch medaka_annotated.vcf.gz
      touch medaka_annotated.vcf.gz.tbi
    else

      cat genome.fasta | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  filename=("reference_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'

      MEDAKA_MODEL=`echo ${params.medaka_model} | sed s'/_variant//'g` # remove _variant form base medaka model
      WSZYSTKIE_CHR=()
      # LONG_CHR lista segmentow ktorych dlugosc wynosi ponad 1000
      # Dla nich wartosc --chunk_len i --chunk_ovlp jest mnozona * 2.5
      LONG_CHR=()
      ALL_VCF=()

      for PLIK in `ls ref*fasta`; do
        SEGMENT_NAME=`basename \${PLIK}  ".fasta" | cut -d "_" -f2-`
        WSZYSTKIE_CHR+=("\${SEGMENT_NAME}")
        LENGTH=`cat reference_\${SEGMENT_NAME}.fasta | grep -v ">" | fold -w1 | wc -l`
        if [ \${LENGTH} -gt 1000 ]; then
          LONG_CHR+=("\${SEGMENT_NAME}")
        fi
      done

      for SEGMENT in \${WSZYSTKIE_CHR[@]}; do
        SEGMENT_IN_BAM=`head -1 reference_\${SEGMENT}.fasta | tr -d ">"` # SEGMENT NAME IN BAM FILE MIGHT CONTAIN EXTRA SYMBOLS LIKE "." or "\"
        # Extract from bam reads assiciated with  specific segment
        samtools view -o initial_\${SEGMENT}.bam trimmed.bam "\${SEGMENT_IN_BAM}"
        samtools index initial_\${SEGMENT}.bam
       
        if [[ \${LONG_CHR[@]} =~ \${SEGMENT} ]]; then
          # This values were used by influenza not sure how SARS or RSV will react
          CHUNK_OVERLAP=`echo "${params.medaka_chunk_overlap}" | awk '{print int(2.5 * \$0)}' ` 
          CHUNK_LEN=`echo "${params.medaka_chunk_len}" |  awk '{print int(2.5 * \$0)}' `
        else
          CHUNK_OVERLAP=`echo "${params.medaka_chunk_overlap}" | awk '{print int(1 * \$0) }'`
          CHUNK_LEN=`echo "${params.medaka_chunk_len}" | awk '{print int(1 * \$0) }'`
        fi

        medaka inference --model \${MEDAKA_MODEL} \
                        --threads ${params.threads} \
                        --chunk_len \${CHUNK_LEN} \
                        --chunk_ovlp \${CHUNK_OVERLAP} \
                        initial_\${SEGMENT}.bam \${SEGMENT}.hdf

        medaka vcf \${SEGMENT}.hdf reference_\${SEGMENT}.fasta medaka_\${SEGMENT}.vcf
        bgzip medaka_\${SEGMENT}.vcf; tabix medaka_\${SEGMENT}.vcf.gz
        ALL_VCF+=("medaka_\${SEGMENT}.vcf.gz")

      done
    
      # merging sub vcf-s
      bcftools concat -a "\${ALL_VCF[@]}" | bcftools sort >> medaka_initial.vcf
    
      medaka tools annotate medaka_initial.vcf genome.fasta pre_trimmed.bam medaka_annotated.vcf
      bgzip medaka_annotated.vcf; tabix medaka_annotated.vcf.gz

      QUAL=`echo ${params.second_round_pval} | awk '{print int(10*-log(\$1)/log(10))}'` 

      bcftools filter -O z -o medaka_annotated_filtered.vcf.gz -i "GQ > \${QUAL} && DP >= ${params.min_cov}" medaka_annotated.vcf.gz
      tabix medaka_annotated_filtered.vcf.gz

    fi
    """
}
