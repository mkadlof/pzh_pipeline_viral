process filtering {
    tag "filtering:${sampleId}"
    container  = params.main_image
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_with_index), val(QC_status), path(primers), path(pairs)
    output:
    tuple val(sampleId), path('to_clip_sorted.bam'), path('to_clip_sorted.bam.bai'), path(primers), path(pairs), val(QC_status), emit: one_amplicon_primers_and_QC
    tuple val(sampleId), path('Primer_usage.txt'), emit: json
    
    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    # Custom filtering for influenza
    # The script takes in sequence:
    # - BAM for filtering and downsampling
    # - Primer scheme
    # - Target coverage per segment (integer)
    # - Minimum read mapping quality (integer)
    # - Reference genome sequences in FASTA (all segments)

    if [ ${QC_status} == "nie" ]; then
      touch to_clip_sorted.bam
      touch to_clip_sorted.bam.bai
      touch Primer_usage.txt
    else
      simple_filter_illumina_INFL.py ${bam} ${primers} ${params.max_depth} ${params.min_mapq} ${params.length} ${ref_genome_with_index[final_index]}
    fi
    """
}

process filtering_nanopore {
    // Nanopore filtering and subsequent primer masking does not require ivar specific pairs file
    tag "filtering:${sampleId}"
    container  = params.main_image
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), path(bai),  path("ref_genome.fasta"), path("primers.bed"), val(QC_status)

    output:
    tuple val(sampleId), path('to_classical_masking.bam'), path('to_classical_masking.bam.bai'), path("ref_genome.fasta"), path("primers.bed"), val(QC_status), emit: to_normal_masking
    tuple val(sampleId), path('to_classical_masking.bam'), path('to_classical_masking.bam.bai'), emit: only_bam
    tuple val(sampleId), path('Primer_usage.txt'), emit: json
    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch to_classical_masking.bam
      touch to_classical_masking.bam.bai
      touch Primer_usage.txt
      QC_exit="nie"
    else
      MASKING_CAP=`echo "${params.mask} + 10" | bc -l` #  do jakiej maksymalnej warotosci podbijac coverage w regionach w ktorych brakowalo oczekiwanych jedno-amplikonowych odczytow
      ALIGNMENT_LENGTH_MAX=1.2 # Ten parametr mowi ze odczyt moze miec 120% dlugosci jego alignmentu z sekwencja referencna
      ALIGNMENT_LENGTH_MIN=0.49 # Ten parametr mowi ze dlugosc regionu mapujacego sie na genom musi stanowic co najmniej 49% dlugosci odczytu 
      # Oba parametry pozwalaja nam odfiltrowac odczyty niepewne na podstawie jaka ich czesc mapuje sie na referencje
      # Odrzucamy zarowno odczyty bardzo dlugie, ktore mapuja sie tylko fragmentarczynie na referencje
      # Jak i takie ktore niezaleznie od swoeje dlugosci maja slabe mapowania 

      python /home/bin/infl/simple_filter_nanopore_INFL_ekstralayer_EQA2024.py ${bam} primers.bed  ${params.bed_offset} ${params.max_depth} ${params.min_mapq} ${params.length} \${ALIGNMENT_LENGTH_MAX} \${ALIGNMENT_LENGTH_MIN} ref_genome.fasta ${params.window_size}

      # W grypie nie ma readow overshot, mergow kilku amplikonow powstaje jeden plik ktory idzie do maskowania
      # Jest on juz tez posortowany w moim skrypcie
      mv to_clip_sorted.bam to_classical_masking.bam 
      samtools index to_classical_masking.bam
    fi
    """
}
