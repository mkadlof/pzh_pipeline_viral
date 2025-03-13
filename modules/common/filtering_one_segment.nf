process filtering {
    tag "filtering:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(primers), path(pairs)

    output:
    tuple val(sampleId), path('first_pass_sorted.bam'), path('first_pass_sorted.bam.bai'), path(primers), path(pairs), val(QC_status), emit: one_amplicon_primers_and_QC
    tuple val(sampleId), path('two_amplicons_sorted.bam'), emit: two_amplicon_only
    tuple val(sampleId), path('Primer_usage.txt'), emit: json
    // tuple val(sampleId), path('Statystyki.txt')

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch first_pass_sorted.bam
      touch first_pass_sorted.bam.bai
      touch two_amplicons_sorted.bam
      touch Primer_usage.txt
      QC_exit="nie"
    else
      bed_offset=0
      length=`echo "${params.length} - 20" | bc -l`
      equal_depth=`echo "${params.max_depth} / 2" | bc -l | awk '{print int(\$0)}'` 
      simple_filter_illumina_one_segment.py ${bam} ${primers} \${bed_offset} \${length} ${params.min_mapq} ${params.window_size} \${equal_depth}
      samtools index first_pass_sorted.bam
    
      if [ -e two_amplicons_sorted.bam ]; then
          cp two_amplicons_sorted.bam tmp.bam
      else
          samtools view -b -o two_amplicons_sorted.bam -L /dev/null first_pass_sorted.bam
      fi
      # Komentarz do json-a - Primer_usage.txt zawiera plik z 3 kolumanmi - nazwa segmentu, nazwa primera, uzyciem
    fi
    """
}

process filtering_nanopore {
    // Nanopore filtering and subsequent primer masking does not require ivar specific pairs file
    tag "filtering:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    // emit bam_and_genome_and_primers z minimap2 
    tuple val(sampleId), path(bam), path(bai),  path("ref_genome.fasta"), path("primers.bed"), val(QC_status)

    output:
    tuple val(sampleId), path('to_classical_masking.bam'), path('to_classical_masking.bam.bai'), path("ref_genome.fasta"), path("primers.bed"), val(QC_status), emit: to_normal_masking
    tuple val(sampleId), path('to_overshot_masking.bam'), path('to_overshot_masking.bam.bai'), path("ref_genome.fasta"), path("primers.bed"), val(QC_status), emit: to_overshot_masking
    tuple val(sampleId), path('to_classical_masking.bam'), path('to_classical_masking.bam.bai'), emit: only_bam
    tuple val(sampleId), path('Primer_usage.txt'), emit: json
    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch to_classical_masking.bam
      touch to_classical_masking.bam.bai
      touch to_overshot_masking.bam
      touch to_overshot_masking.bam.bai
      touch Primer_usage.txt
    else
      MASKING_CAP=`echo "${params.mask} + 10" | bc -l` #  do jakiej maksymalnej warotosci podbijac coverage w regionach w ktorych brakowalo oczekiwanych jedno-amplikonowych odczytow
      # jako ze korzystamy z odczytow o niejasnym pochodzeniu (najczesciej odczty z fuzji amplikonow)
      # to dobijamy tylk odo wartosci tak aby podbic zeby region nie byl maskowany
      equal_depth=`echo "${params.max_depth} / 2" | bc -l | awk '{print int(\$0)}'`  
      simple_filter_nanopore_final_with_windowstep.py ${bam} primers.bed ${params.bed_offset} \${equal_depth} ${params.length} ${params.min_mapq} ${params.extra_bed_offset} \${MASKING_CAP} ${params.window_size}

      # Skrypt wyzej zwraca bardzo duzo plikow, niestety aktualnie ich powstawanie jest zalezne od danych 
      # (zawsze bedzie reads_inner_strict.bam, reszta jest opcjonalna)

      TO_MERGE_CLASSICAL=()
      TO_MERGE_EXTRA=()
      
      # Odczyty najlepsze
      if [ -e reads_inner_strict.bam ]; then
        TO_MERGE_CLASSICAL+=('reads_inner_strict.bam')
      else
        # To by oznaczalo ze cos jest bardzo zle
        touch dummy.txt
      fi

      # Odczyty mapujace sie na jeden primer, moga byc ale nie zawsze
      if [ -e reads_partial_strict_sort.bam ]; then
        TO_MERGE_CLASSICAL+=('reads_partial_strict_sort.bam')
      fi

      # Odczyty tylko w protokole midnight, truden do okreselenie pochodzenie
      # Uzywane do podbicia pokrycia w sytuacji ratunkowej
      if [ -e reads_smieci_sorted.bam ]; then
        TO_MERGE_CLASSICAL+=('reads_smieci_sorted.bam')
      fi 
      
      # Laczenie plikow ktore beda maskowane z ta sama tolerancja, sortowanie i indeksowanie 
      samtools merge -o to_classical_masking_initial.bam \${TO_MERGE_CLASSICAL[@]}
      samtools sort -o to_classical_masking.bam to_classical_masking_initial.bam
      samtools index to_classical_masking.bam

      # Odczyty pochodzace z dwoch amplikonow sa w oddzielnych plikach nalezy je polaczyc w jeden 
      # plik razem z plikiem reads_inner_strict.bam zawierajacym "dobre" odczyty
      
      if [ `ls -l  reads_two_amplicons*bam | wc -l` -gt 0 ];then
         for PLIK in `ls reads_two_amplicons*bam`; do
           TO_MERGE_EXTRA+=("\${PLIK}")
         done 
      fi

      if [ -e reads_overshot.bam ]; then
        TO_MERGE_EXTRA+=("reads_overshot.bam")
      fi

      
      #  Laczenie plikow ktore beda maskowane z wyzsza tolerancja, ich sortowanie i indeksowanie
      # Tych plikow moze nie byc stas sa pod if-em
      if [ "\${#TO_MERGE_EXTRA[@]}" -gt 0 ]; then
        samtools merge -o to_overshot_masking_initial.bam \${TO_MERGE_EXTRA[@]}
        samtools sort -o to_overshot_masking.bam to_overshot_masking_initial.bam
        samtools index to_overshot_masking.bam
      else
        # Stworzmy pliki bam 
        samtools view -H -o to_overshot_masking.bam reads_inner_strict.bam
        samtools index to_overshot_masking.bam
      fi

    fi
    """
}
