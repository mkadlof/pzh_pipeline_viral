process reassortment {
    tag "reassortment:${sampleId}"
    container  = params.main_image 
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "hybrid_genome.fasta"
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path('subtype_mean_coverage_each_segment.txt'), path('subtype_scores_each_segment.txt'), val(REF_GENOME_ID_entry), val(QC_status)
    // path(genomes)
    // path(primers)

    output:
    tuple val(sampleId), path("hybrid_genome.fasta*"), path("hybrid_primers.bed"), path("pairs.tsv"), env(REF_GENOME_ID), env(QC_exit), emit: all
    tuple val(sampleId), path("hybrid_genome.fasta"), path("hybrid_primers.bed"), env(QC_exit), emit: all_nanopore
    tuple val(sampleId), path("hybrid_genome.fasta*"), path("hybrid_primers.bed"), env(QC_exit), emit: to_bwa
    tuple val(sampleId), path("hybrid_primers.bed"), path("pairs.tsv"), emit: primers_and_pairs
    tuple val(sampleId), path("hybrid_primers.bed"), emit: primers
    tuple val(sampleId), path("hybrid_genome.fasta"), emit: only_genome // indelqual module requires a variable not a tupple
    tuple val(sampleId), path("reassortment.json"), emit: json
    tuple val(sampleId), path("hybrid_genome.fasta"), path("genes.gtf"), emit: to_snpeff
    script:
    """
    find_segement_position() {
        local VALUES=("\${!1}")
        local hit=\$2
        local i=1
        for ela in \${VALUES[@]}
        do
            if [ "\$ela" == "\$hit" ]; then
                echo \${i}
                return 0
            else
                ((i++))
            fi
        done
    }
    # REF_GENOME_ID is defined here again, 'cause I am lazy
    if [ ${QC_status} == "nie" ]; then
      touch hybrid_primers.bed
      touch hybrid_genome.fasta
      touch hybrid_genome.fasta.fai
      touch pairs.tsv
      REF_GENOME_ID="unk"
      QC_exit="nie"
      touch genes.gtf
      # json section 
      ERR_MSG="This module recieved failed QC status"
      STATUS="nie"
      influenza_reassortment_parser.py --status "\${STATUS}" \
                                        --error "\${ERR_MSG}" \
                                        --output reassortment.json

    else 
      REF_GENOME_ID="${REF_GENOME_ID_entry}"

      ALL_GENOMES=(`ls /home/data/infl/genomes`)
      ALL_SEGMENTS=(PB2 PB1 PA HA NP NA MP NS)
      genomes="/home/data/infl/genomes" # path to genomes WITHIN container
      primers="/home/data/infl/primers" # path to primers WITHIN container

      cat \${genomes}/\${REF_GENOME_ID}/\${REF_GENOME_ID}.fasta | \
      awk -v ID=\${REF_GENOME_ID} '{if (substr(\$0, 1, 1)==">") {filename=(ID"_"substr(\$0,2) ".fasta"); print \$0"_"ID >> filename } else {print toupper(\$0)  >> filename}}'

      # For each segment, we analyze the file subtype_scores_each_segment.txt. This file holds a matrix where
      # rows represent subtypes, and columns mapping scores to segments from a given subtype .

      SEGMENTS_REASSORTMENT=()
      FOUND_SUBTYPES_REASSORTMENT=()
      FOUND_SUBTYPES_SCORE_RATIO=()
      FOUND_SUBTYPES_SEQ_SIMILARITY=()
      FOUND_SUBTYPES_COUNTS=()

      # These 5 lists should hold, in order: the segment identifier, the subtype with the highest
      # alignment score  the score of the expected subtype / score of the best subtype, the
      # similarity between the expected subtype and the found subtype, and the average coverage on
      # the best segment. These list are only used for internal QC and are printed to intermediate.txt file
   
      # Basically "expected" subtype is subtype identified based on mapping score for HA and NA only
      # "Best" is the subtype with highest mapping score for that segment  


      for segment in \${ALL_SEGMENTS[@]}; do
        SEGMENTS_REASSORTMENT+=(\${segment})

        # We need to extract data from a file so first we determine in which column the data
        # for a given segment is located.
        LIST_OF_IDS=`cat subtype_scores_each_segment.txt | head -1`
        SEGMENT_POS=`find_segement_position LIST_OF_IDS[@] \${segment}`

        # Note we moved from simple counting to mean coverage at some point
        # So the names of variables can be misleading, we should fix that later
        LIST_OF_IDS_COUNTS=`cat subtype_mean_coverage_each_segment.txt | head -1`
        SEGMENT_POS_COUNTS=`find_segement_position LIST_OF_IDS_COUNTS[@] \${segment}`

        # name of the subtype with the best result for this segment
        SEGMENT_best=`cat subtype_scores_each_segment.txt | sort -rnk\${SEGMENT_POS} | head -1 | cut -f1`

        # what is the alignment score for this segment from this subtype
        SEGMENT_best_score=`cat subtype_scores_each_segment.txt | sort -rnk\${SEGMENT_POS} | head -1 | cut -f\${SEGMENT_POS}`

        # what is the score for this segment if taking a segemnt from  mapped HA/NA?
        SEGMENT_expected_score=`cat subtype_scores_each_segment.txt | grep \${REF_GENOME_ID} | cut -f \${SEGMENT_POS}`

        SEGMENT_best_counts=`cat subtype_mean_coverage_each_segment.txt | sort -rnk\${SEGMENT_POS_COUNTS} | head -1 | cut -f\${SEGMENT_POS_COUNTS}`
        FOUND_SUBTYPES_COUNTS+=(\${SEGMENT_best_counts})

        if [ \${SEGMENT_best} == \${REF_GENOME_ID} ]; then
            # I am not doing anything, subtype which segment obtained highest mapping scor
            # is identical to subtype determined using HA/NA segments only
            FOUND_SUBTYPES_REASSORTMENT+=(\${REF_GENOME_ID})
            FOUND_SUBTYPES_SCORE_RATIO+=(1)
            FOUND_SUBTYPES_SEQ_SIMILARITY+=(1)

            # Creating a "hybrid" for the genome and BED with primers. We do this N times because
            # if any other segment is not a reassortment, the FASTA files will "grow" with copies
            # of the same sequence.

            cat \${genomes}/\${REF_GENOME_ID}/\${REF_GENOME_ID}.fasta | \
                awk -v ID=\${REF_GENOME_ID} '{if (substr(\$0, 1, 1)==">") {filename=("regular_"ID"_"substr(\$0,2) ".fasta");  print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
            
            cat regular_\${REF_GENOME_ID}_chr?_\${segment}.fasta >> hybrid_genome.fasta
            
            cat \${primers}/\${REF_GENOME_ID}/\${REF_GENOME_ID}_primers.bed | \
                grep \${segment} >> hybrid_primers.bed

           cat /opt/snpEff/data/\${REF_GENOME_ID}/genes.gtf | \
                grep \${segment} >> genes.gtf
 
            # remove intermediate files
            rm regular_\${REF_GENOME_ID}_chr*
        else
            # subtype for which analyzed segment had the highest mapping score
            # is different from subtype determined based on the HA/NA combination

            FOUND_SUBTYPES_REASSORTMENT+=(\${SEGMENT_best})
            RATIO=`echo "\${SEGMENT_expected_score} / \${SEGMENT_best_score}" | bc -l`
            FOUND_SUBTYPES_SCORE_RATIO+=(\$RATIO)
            
            # To determine if we observe "true" reassortment we need to check
            # what is the similarity between segments form "best" and "expected" subtypes

            cat \${genomes}/\${SEGMENT_best}/\${SEGMENT_best}.fasta | \
                awk -v ID=\${SEGMENT_best} '{if (substr(\$0, 1, 1)==">") {filename=(ID"_"substr(\$0,2) ".fasta");  print \$0"_"ID >> filename } else {print toupper(\$0) >> filename}}'

            # And also a "regular" FASTA for the hybrid
            cat \${genomes}/\${SEGMENT_best}/\${SEGMENT_best}.fasta | \
                awk -v ID=\${SEGMENT_best} '{if (substr(\$0, 1, 1)==">") {filename=("regular_"ID"_"substr(\$0,2) ".fasta");  print \$0 >> filename } else {print toupper(\$0) >> filename}}'

            cat \${REF_GENOME_ID}_chr?_\${segment}.fasta \${SEGMENT_best}_chr?_\${segment}.fasta >> tmp.fa
            
            # Run NW and extract score
            SEGMENT_alignment_score=`run_nw.py tmp.fa | bc -l`
            FOUND_SUBTYPES_SEQ_SIMILARITY+=(\${SEGMENT_alignment_score})

            # Empirically true reassortment  is when we analyzed segments have diverse sequence (determined by NS)
            # the mapping score ratio of segement form expected subtype is less than 0.6 of score achieved by segment
            # from "best" subtype and avarage coverage predicted for segment from "best" subtype is at least 150% 
            # of minimum coverage
            
            if awk "BEGIN {exit !(\${SEGMENT_alignment_score} < 0.9 && \${RATIO} < 0.6 && \${SEGMENT_best_counts} >= (${params.min_cov} * 1.1))}"; then
              # echo "Reassortment detected for the segment \${segment}: \${REF_GENOME_ID} -> \${SEGMENT_best}\tAverage coverage is \${SEGMENT_best_counts}"

              cat regular_\${SEGMENT_best}_chr?_\${segment}.fasta >> hybrid_genome.fasta

              cat /home/data/infl/primers/\${SEGMENT_best}/\${SEGMENT_best}_primers.bed | \
                  grep \${segment} >>  hybrid_primers.bed

              cat /opt/snpEff/data/\${SEGMENT_best}/genes.gtf | \
                grep \${segment} >> genes.gtf

            else
                # Probbably not a true reassortment we use segment from the "expected" subtype
                cat \${genomes}/\${REF_GENOME_ID}/\${REF_GENOME_ID}.fasta | \
                    awk -v ID=\${REF_GENOME_ID} '{if (substr(\$0, 1, 1)==">") {filename=("regular_"ID"_"substr(\$0,2) ".fasta");  print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
                cat regular_\${REF_GENOME_ID}_chr?_\${segment}.fasta >> hybrid_genome.fasta

                cat \${primers}/\${REF_GENOME_ID}/\${REF_GENOME_ID}_primers.bed | \
                     grep \${segment} >> hybrid_primers.bed

                cat /opt/snpEff/data/\${REF_GENOME_ID}/genes.gtf | \
                grep \${segment} >> genes.gtf

                rm regular_\${REF_GENOME_ID}_*
            fi
            rm tmp.fa
            rm *\${SEGMENT_best}_chr*.fasta
        fi
      done

      # intermediate file for QC checks
      rm \${REF_GENOME_ID}_chr*.fasta
      echo \${SEGMENTS_REASSORTMENT[@]} >> intermediate.txt
      echo \${FOUND_SUBTYPES_REASSORTMENT[@]} >> intermediate.txt
      echo \${FOUND_SUBTYPES_SCORE_RATIO[@]} >> intermediate.txt
      echo \${FOUND_SUBTYPES_SEQ_SIMILARITY[@]} >> intermediate.txt
      echo \${FOUND_SUBTYPES_COUNTS[@]} >> intermediate.txt

      REFERENCE_GENOME_FASTA="hybrid_genome.fasta"
      bwa index \${REFERENCE_GENOME_FASTA}
      PRIMERS="hybrid_primers.bed"
      # Pairs file is subtype independent for influenza
      cp /home/data/infl/primers/pairs.tsv .
      QC_exit="tak"

      #json section
      influenza_reassortment_parser.py --status "tak" \
                                        --output reassortment.json \
                                        --input_file intermediate.txt \
                                        --subtype \${REF_GENOME_ID} \
                                        --mapping /home/data/infl/strains_name.txt
   fi #koniec if-a na QC
   """
}
