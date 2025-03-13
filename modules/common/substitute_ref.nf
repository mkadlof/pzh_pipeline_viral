process substitute_ref_genome {
    // This process is used to substitute reference genome from second stept of genome generation
    // in nanopore to the originally used reference genome and to integrate nanopore and illumina I/O
    // It also serves as a QC switch to mimic functionality of introduce_SV_with_manta module from illumina
    tag "genome_substitution:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status), path('original.fasta')
    output:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path('original_reference.*'), env(QC_status_exit), emit: fasta_refgenome_and_qc
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), env(QC_status_exit), emit: fasta_and_qc // for novel snpeff
    script:
   """
   if [ ${QC_status} == "nie" ]; then
     touch output_consensus_masked_SV.fa
     touch original_reference.fasta
     QC_status_exit="nie"
   else
     cp original.fasta original_reference.fasta
     bwa index original_reference.fasta
    
     NUMBER_OF_N=`cat output_consensus_masked_SV.fa  | grep -v ">" | fold -w1 | sort | uniq -c | grep N | cut -d " " -f3`
     SEQ_LENGTH=`cat output_consensus_masked_SV.fa | grep -v ">" | fold -w1 | wc -l`
     if [ -z \${NUMBER_OF_N} ]; then
       # No Ns in a sequence
       QC_status_exit="tak"
     else
       if [ `awk -v n="\${NUMBER_OF_N}" -v total="\${SEQ_LENGTH}" 'BEGIN {wynik=n/total; if (wynik < 0.9) print "1"; else print "0"}'` -eq 1 ]; then
         # less than 90% of Ns in sequence we can continue ...
         QC_status_exit="tak"
       else
         # ore than 90% of Ns
         QC_status_exit="nie"
       fi
     fi
   fi
   """
}
