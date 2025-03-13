process picard_downsample {
    tag "picard:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path('downsample.bam'), path('downsample.bam.bai'), path(ref_genome_with_index), env(QC_exit), emit: to_manta

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch downsample.bam
      touch downsample.bam.bai
      QC_exit="nie"
    else

      ILE=`samtools view ${bam} | wc -l `
      ILE=`echo "\${ILE} + 2" | bc -l`
      NORM=`echo "${params.max_number_for_SV}/\${ILE}" | bc -l `
      if [ `awk -v norm="\${NORM}" 'BEGIN {if(norm>0.99) {print 1} else {print 0}}'` -eq 1 ]; then 
        NORM=0.99
      fi
      echo \${NORM}
      suffix=`echo "${params.max_number_for_SV}/1000" | bc -l | awk '{print int(\$0)}'`

      java -jar /opt/picard/picard.jar PositionBasedDownsampleSam \
                                      --INPUT ${bam} \
                                      --OUTPUT downsample.bam -F \${NORM}
      samtools index downsample.bam
      NO_READS=`samtools view downsample.bam | wc -l`
      if [ "\${NO_READS}" -lt 10000 ]; then
        QC_exit="nie"
      else
        QC_exit="tak"
      fi # koniec if na za malo odczytow w downsampled bam
    
    fi # koniec if na initial QC_status
    """
}

process picard_downsample_multisegment {
    tag "picard:${sampleId}"
    // Generalizacja moduly dla picarda dla jednego semgnetu
    // QC_exit jest "nie" jesli WSZYSTKIE segmenty zwroca bam w ktorym MEDIANA pokrycie jest ponizej 50 
    // W ten sposob unikamy probek w ktorym calosc odczytow to absurdalny spike pokrycia w waskim regionie
    // segmentu przy "pustych" czesciach resty segmentu
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_with_index), val(QC_status)

    output:
    tuple val(sampleId), path('downsample*.bam'), path('downsample*.bam.bai'), path(ref_genome_with_index), path("mediana_per_segment.txt"), env(QC_exit), emit: to_manta

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch downsample.bam
      touch downsample.bam.bai
      touch mediana_per_segment.txt
      QC_exit="nie"
    else
      MEDIAN_ALL=()
      # Spli initial bam into per segments bam
      SEGMENTS_NAME=(`samtools view ${bam} | cut -f3 | sort | uniq | tr "\n" " "`)
      for segment in \${SEGMENTS_NAME[@]}; do
        # split multisegment bam into segment unique bams
        segment_clean=`echo \${segment} | tr "." "_" | tr "/" "_"` # segment clean 3ma nazwe w ktorej symbole "." (genom sars-a), "/" (genom RSV) sa zamienione na "_" w nazwch plikow,
        # docelowo zawartosc takiego pliku 3ma oryginalna nazwe segmentu

        samtools view -o initial_\${segment_clean}.bam ${bam} "\${segment}"

        ILE=`samtools view  initial_\${segment_clean}.bam | wc -l `
        ILE=`echo "\${ILE} + 2" | bc -l`
        NORM=`echo "${params.max_number_for_SV}/\${ILE}" | bc -l `

        if [ `awk -v norm="\${NORM}" 'BEGIN {if(norm>0.99) {print 1} else {print 0}}'` -eq 1 ]; then
          NORM=0.99
        fi

        java -jar /opt/picard/picard.jar PositionBasedDownsampleSam \
                                         --INPUT initial_\${segment_clean}.bam \
                                         --OUTPUT downsample_\${segment_clean}.bam -F \${NORM}
        samtools index downsample_\${segment_clean}.bam
     
        # Liuczymy srednie pokrycie w seg 
        bedtools genomecov -d -ibam downsample_\${segment_clean}.bam > downsample_\${segment_clean}.bed
        MEDIANA_SEGMENT=`sort -nk3 downsample_\${segment_clean}.bed | awk '{a[i++]=\$3} END {print a[int(i/2)]}'`
        echo -e "\${segment_clean}\t\${MEDIANA_SEGMENT}" >> mediana_per_segment.txt
        MEDIAN_ALL+=(\${MEDIANA_SEGMENT})
      done
      # Now check the medians in all the segments if at least one of the values is above threshols 
      # change QC_exit to "tak"
      QC_exit="nie"
      for val in \${MEDIAN_ALL[@]}; do
        # UWAGA MEDIANA W SEGMENCIE JEST POWYZEJ 50 WIEC PROBUJEMY CALOWAC SV !!!
        if [ \${val} -ge 50 ]; then
          QC_exit="tak"
        fi
      done


    fi # koniec if na initial QC_status
    """
}
