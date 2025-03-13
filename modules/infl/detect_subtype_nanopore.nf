process detect_subtype_nanopore {
    // This process is virtually identical to illumina
    // but we use minimap2 instead of bwa and reads variable is not a tuple
    tag "detect_subtype:${sampleId}"
    container  = params.main_image
    cpus params.threads
    memory "40 GB"
    input:
    tuple val(sampleId), path(reads), val(QC_STATUS)
    output:
    tuple val(sampleId), path("subtype_mean_coverage_each_segment.txt"), path("subtype_scores_each_segment.txt"), env(REF_GENOME_ID), val(QC_STATUS), emit: segments_scores
    tuple val(sampleId), env(REF_GENOME_ID_MINI), emit: subtype_id
    tuple val(sampleId), path(reads), env(REF_GENOME_ID_MINI), val(QC_STATUS), emit: to_freyja

    script:
    """
    # Functions
    run_minimap() {
    for GENOME in "\${@}"; do
        minimap2 -a -x map-ont -t ${task.cpus} -o \${GENOME}.sam /home/data/infl/genomes/\${GENOME}/\${GENOME}.fasta ${reads}
        samtools view -bS -F 2052 \${GENOME}.sam > \${GENOME}.bam
        samtools sort -o \${GENOME}_sorted.bam \${GENOME}.bam
        mv \${GENOME}_sorted.bam \${GENOME}.bam
        samtools index \${GENOME}.bam
        rm \${GENOME}.sam
    done

   }

    find_max() {
        local GENOMES=("\${!1}")
        local VALUES=("\${!2}")
        local max=0
        local max_index=0
        local i=0
        for value in \${VALUES[@]}
        do
            echo \${GENOMES[\${i}]} \${value} >> subtype_scores.txt
            if [ \${value} -gt \${max} ]; then
                max=\${value}
                max_index=\${i}
            fi
            ((i++))
        done
        echo \${GENOMES[\${max_index}]}
    }

    # Main

    if [ ${QC_STATUS} == "nie" ]; then
      touch subtype_mean_coverage_each_segment.txt
      touch subtype_scores_each_segment.txt
      REF_GENOME_ID="unk"
      REF_GENOME_ID_MINI="unk"
    else

      #KNOWN_VARIANTS='H1N1 H3N2 H4N6 H5N2 H5N1 H5N6 H5N8 H6N1 H7N9 H9N2 Yamagata Victoria UNK'
      ALL_GENOMES=(`ls /home/data/infl/genomes`)
      ALL_SEGMENTS=(PB2 PB1 PA HA NP NA MP NS)

      echo -e "id \${ALL_SEGMENTS[@]}" | tr " " "\t" >> subtype_mean_coverage_each_segment.txt
      echo -e "id \${ALL_SEGMENTS[@]}" | tr " " "\t" >> subtype_scores_each_segment.txt

      if [ ${params.variant} == 'UNK' ] ;then
        # User does not known the varaint of his influenza sample
        # we need to predict it ourselves
        # We map our sequences to each of the references and count the number of reads mapping to
        # the HA and NA segments. The reference with most reads is selceted.

        run_minimap \${ALL_GENOMES[@]}

        for GENOME in \${ALL_GENOMES[@]}; do
            alignment_score+=(`get_alignment_score.py \${GENOME}.bam chr6_NA,chr4_HA`)
            echo -e "\${GENOME}\t`get_alignment_score_all_segments.py \${GENOME}.bam`" >> subtype_scores_each_segment.txt
            MEAN_COVERAGE=()
            for SEGMENT in \${ALL_SEGMENTS[@]}; do
                MEAN_COVERAGE+=(`bedtools genomecov -ibam \${GENOME}.bam  -bga | \
                                 grep -w chr[0-9]_\${SEGMENT} | \
                                 awk 'BEGIN {segment_length=0; read_count=0} {segment_length+=(\$3-\$2); read_count+=\$4 * (\$3-\$2)} END {print read_count/segment_length}'`)
            done
            echo -e "\${GENOME} \${MEAN_COVERAGE[@]}" | tr " " "\t" >> subtype_mean_coverage_each_segment.txt
        done

        result=`find_max ALL_GENOMES[@] alignment_score[@]`
        # echo \${result}
        result_mini=`echo \${result} | cut -d "_" -f1`

      else
        # The procedure is mostly identical to "UNK", however we only analyze genomes
        # for a particulat user-specified subtype, in our case for some, popular subtypes
        # we have references for main clades
        SAMPLE_GENOMES=()

        for genome in "\${ALL_GENOMES[@]}"; do
            if [[ "\${genome}" == *"\${variant}"* ]]; then
                    SAMPLE_GENOMES+=(\${genome})
            fi
        done

        run_minimap \${SAMPLE_GENOMES[@]}

        alignment_score=()
        for GENOME in \${SAMPLE_GENOMES[@]}; do
            alignment_score+=(`get_alignment_score.py \${GENOME}.bam chr6_NA,chr4_HA`)
            echo -e "\$GENOME\t`get_alignment_score_all_segments.py \${GENOME}.bam`" >> subtype_scores_each_segment.txt
            MEAN_COVERAGE=()
            for SEGMENT in \${ALL_SEGMENTS[@]}; do
                MEAN_COVERAGE+=(`bedtools genomecov -ibam  \${GENOME}.bam  -bga | grep -w chr[0-9]_\${SEGMENT} | awk 'BEGIN {segment_length=0; read_count=0} {segment_length+=(\$3-\$2); read_count+=\$4 * (\$3-\$2)} END {print read_count/segment_length}'`)
            done
            echo -e "\${GENOME} \${MEAN_COVERAGE[@]}" | tr " " "\t" >> subtype_mean_coverage_each_segment.txt
        done

        result=`find_max SAMPLE_GENOMES[@] alignment_score[@]`
        result_mini=`echo \${result} | cut -d "_" -f1`

      fi # koniec if-a na nieznany wariant
      REF_GENOME_ID=\${result}
      REF_GENOME_ID_MINI=`echo \${result} | cut -d "_" -f1`
      rm *bam*
    fi # koniec if-a na zle QC
    """
}
