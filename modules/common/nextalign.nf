process nextalign {
    tag "nextalign:${sampleId}"
    container = params.main_image
    cpus 1
    memory "20 GB"
    // publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "*.fasta"

    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), val(QC_status), val(SAMPLE_SUBTYPE)

    output:
    // Will depend on organism for infl this module produce output for antiviral resistance and modeller
    // For SARS2 and RSC it produce only output to modeller
    tuple val(sampleId), path("nextalign_gene*fasta"), val(QC_status), emit: to_modeller
    tuple val(sampleId), path("to_resistance.fasta"), path("sample_M2.fasta"), val(QC_status), val(SAMPLE_SUBTYPE),  optional: true, emit: to_resistance

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch to_resistance.fasta
      touch sample_M2.fasta
      touch nextalign_gene_dummy.fasta
    else
      # all data for nextalign are predefined and sotred in /home/data/common/nextalign
      
     if [ ${params.species} == "SARS-CoV-2" ]; then

        # -r reference genomic sequence; -m gff file with protein location; 
        # -g protein name '.' for all proteins in gff;
        #  -O output dir; positional argument sequence of you sample
        nextalign run -r /home/data/common/nextalign/MN908947.3/sequences.fa \
                      -m /home/data/common/nextalign/MN908947.3/genes.gff3 \
                      -g "E,M,N,ORF10,orf1a,orf1b,ORF3a,ORF6,ORF7a,ORF8,S" \
                      -O . \
                      output_consensus_masked_SV.fa
        # remove "-" from protein sequences + change their header to include saple id (after "|")
        # TO DO
      elif [ ${params.species} == "RSV" ]; then
        nextalign run -r /home/data/common/nextalign/hRSV_${SAMPLE_SUBTYPE}/sequences.fa \
                      -m /home/data/common/nextalign/hRSV_${SAMPLE_SUBTYPE}/genes.gff3 \
                      -g "NS1,NS2,N,P,M,SH,G,F,M2,M2,L" \
                      -O . \
                      output_consensus_masked_SV.fa
        # remove "-" from protein sequences + change their header to include saple id (after "|")
        # TO DO
      elif [ ${params.species} == "Influenza" ]; then

        cat output_consensus_masked_SV.fa | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  gsub("_SV", "", new_name);  filename=("sample_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
      
        NEXTALIGN_DB="/home/data/common/nextalign"
      
        # For these subtypes we prepared custom files to get sequences of ALL proteins 
        TO_NEXTALIGN='H1N1 H3N2 H4N6 H5N1 H5N2 H5N5 H5N6 H5N8 H6N1 H7N9 H9N2 Yamagata Victoria'
        if [[ \${TO_NEXTALIGN[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
          for FILE in `ls sample_*.fasta`
          do
            SEGMENT=`basename \${FILE} ".fasta" | cut -d "_" -f3 `
            nextalign run -r \${NEXTALIGN_DB}/${SAMPLE_SUBTYPE}/\${SEGMENT}.fasta \
                          -m \${NEXTALIGN_DB}/${SAMPLE_SUBTYPE}/${SAMPLE_SUBTYPE}.gff \
                          -g \${SEGMENT} \
                          -O . \${FILE}

            if [[ \${SEGMENT} == 'NA' || \${SEGMENT} == 'PA' ]]; then
                cat nextalign_gene_\${SEGMENT}.translation.fasta >> to_resistance.fasta
            fi

            if [ \${SEGMENT} == 'MP' ]; then
                # Te jest glupie ale moj skrypt ma hard coded ta nazwe 
                cp sample_chr7_MP.fasta consensus_MP.fasta
                # I hard kodowny naglowek
                sed -i s"/chr7_MP_SV/MP/"g consensus_MP.fasta
                sed -i s"/chr7_MP/MP/"g consensus_MP.fasta
                prep_M2.py ${SAMPLE_SUBTYPE} sample_M2.fasta \${NEXTALIGN_DB} . .
                nextalign run -r M2.fasta \
                              -m \${NEXTALIGN_DB}/${SAMPLE_SUBTYPE}/${SAMPLE_SUBTYPE}.gff \
                              -g M2 \
                              -O . consensus_M2.fasta
               cp nextalign_gene_M2.translation.fasta sample_M2.fasta
            fi
          done
        fi # koniec if-a na dostepnosc bazy nextalign
      fi # koniec if-a na organizm
    fi # koniec if na QC
    """
}
