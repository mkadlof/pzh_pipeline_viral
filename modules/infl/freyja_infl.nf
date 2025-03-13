process freyja_infl {
    // We need to repeat mapping for RSV as Freyja needs a different genome than we used in main pipeline
    // We also "attach" this module not AFTER bwa/minmap2 but NEXT to it
    // We do not mask primers here
 
    tag "freyja:${sampleId}"
    container  = params.main_image
    cpus params.threads
    memory "20 GB"
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"

    input:
    tuple val(sampleId), path(reads), val(SAMPLE_SUBTYPE), val(QC_status)

    output:
    tuple val(sampleId), path('coinfections_freyja.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch coinfections.tsv
      freyja_status="nie"
      ERR_MSG="This module recieved a failed QC status and was not executed"
      echo -e "{\\"status\\":\\"\${freyja_status}\\",
                \\"error_message\\":\\"\${ERR_MSG}\\"}" >> coinfections_freyja.json
    else
      TO_FREYJA='H1N1 H3N2 H5N1 H5N2 H5N5 H5N6 H5N6 Victoria'
      if [[ \${TO_FREYJA[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
        if [[ ${SAMPLE_SUBTYPE} == *H5* ]]; then
          FREYJA_SUBTYPE="H5Nx"
        else
          FREYJA_SUBTYPE="${SAMPLE_SUBTYPE}"
	fi

         mkdir variants_files depth_files demix_files
         cp /home/external_databases/freyja/\${FREYJA_SUBTYPE}/reference.fasta .
         bwa index reference.fasta

         # We map  reads all the reads to Freyja-required genome 
         if [ ${params.machine} == 'Illumina' ]; then  
           bwa mem -t ${task.cpus} -T 30 reference.fasta ${reads[0]} ${reads[1]} | \
           samtools view -@ ${task.cpus} -Sb -f 3 -F 2048 - | \
           samtools sort -@ ${task.cpus} -o mapped_reads.bam -
           samtools index mapped_reads.bam
         elif [ ${params.machine} == 'Nanopore' ]; then
           minimap2 -a -x map-ont -t ${task.cpus} -o tmp.sam reference.fasta ${reads}
           samtools view -@ ${params.threads} -Sb -F 2052 tmp.sam | \
           samtools sort -@ ${params.threads} -o mapped_reads.bam -
           samtools index mapped_reads.bam
         fi


         # We call variants with freyja
 
         freyja variants mapped_reads.bam --minq ${params.freyja_minq} --variants variants_files/test.variants.tsv --depths depth_files/test.depth --ref reference.fasta
         freyja demix variants_files/test.variants.tsv depth_files/test.depth --output demix_files/test.output --confirmedonly --barcodes  /home/external_databases/freyja/\${FREYJA_SUBTYPE}/barcode.csv

         freyja aggregate demix_files/ --output coinfections.tsv
         freyja_status="tak"
         freyja_lineage_1_name=`cat coinfections.tsv  | cut -f3 | tail -1 | cut -d " " -f1 | sed s"|\${FREYJA_SUBTYPE}-||"g`
         freyja_lineage_2_name=`cat coinfections.tsv  | cut -f3 | tail -1 | cut -d " " -f2 |  sed s"|\${FREYJA_SUBTYPE}-||"g`
         freyja_lineage_1_abundance=`cat coinfections.tsv  | cut -f4 | tail -1 | cut -d " " -f1 | awk '{printf "%.2f", \$0}'`
         freyja_lineage_2_abundance=`cat coinfections.tsv  | cut -f4 | tail -1 | cut -d " " -f2 | awk '{printf "%.2f", \$0}'`

         # In case freyja found a single linage
         if [ -z \${freyja_lineage_2_name} ]; then
           freyja_lineage_2_name="unk"
           freyja_lineage_2_abundance=0
         fi

         # In case both linages have the same name
         if [ \${freyja_lineage_1_name} == \${freyja_lineage_2_name} ]; then
           freyja_lineage_2_name="unk"
           freyja_lineage_2_abundance=0
         fi

         echo -e "{\\"status\\":\\"\${freyja_status}\\",
                   \\"freyja_lineage1_name\\":\\"\${freyja_lineage_1_name}\\",
                   \\"freyja_lineage2_name\\":\\"\${freyja_lineage_2_name}\\",
                   \\"freyja_lineage1_abundance\\":\${freyja_lineage_1_abundance},
                   \\"freyja_lineage2_abundance\\":\${freyja_lineage_2_abundance}}" >> coinfections_freyja.json
       else
         ERR_MSG="Freyja module works only with following influenza subtypes H1N1, H3N2, H5Nx and Victoria"
         touch coinfections.tsv
         freyja_status="nie"
         ERR_MSG="This module recieved a failed QC status and was not executed"
         echo -e "{\\"status\\":\\"\${freyja_status}\\",
                   \\"error_message\\":\\"\${ERR_MSG}\\"}" >> coinfections_freyja.json
       fi

    fi

    """
}
