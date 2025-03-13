process vcf_from_fasta {
    tag "vcf_from_fasta:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    // publishDir "${params.results_dir}/${sampleId}/consensus_vcf", mode: 'copy', pattern: "*vcf*"

    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), val(QC_status), path('sequences.fa'), path('genes.gtf')

    output:
    tuple val(sampleId), path("${sampleId}_final.vcf.gz"), path("${sampleId}_final.vcf.gz.tbi"), val(QC_status), path('sequences.fa'), path('genes.gtf'), emit: vcf

    """
    # Skrypt nie poradzi sobie z fasta z wieloma segmentami, wiec analizujemy genome per segment i na koncu sklejamy vcf-a
    
    
     if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_final.vcf.gz
      touch ${sampleId}_final.vcf.gz.tbi
      touch sequences.fa
      touch genes.gtf
    else

      GENOME_FASTA="sequences.fa"
 
      cat \${GENOME_FASTA} | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  filename=("reference_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
      cat output_consensus_masked_SV.fa | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  gsub("_SV", "", new_name);  filename=("sample_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}' 
      LISTA_VCF=()
      for PLIK in `ls reference_*`; do
        SEGMENT=`basename \${PLIK} ".fasta" | cut -d "_" -f2- `
        N=1 # maksymalna przerwa miedzy mutacjami aby je polaczyc w jeden ciag
        vcf_input="/home/data/sarscov2/vcf_template/vcf_template.vcf.gz" # sciezka do pliku vcf, paczka vcf z pythona wymaga pliku vcf aby na podstawie tego schematy samemu tworzyc plik output-owy
        vcf_output="consensus_\${SEGMENT}.vcf" # nazwa pliku vcf dla sekwencji konsensusowej
        sed -i s'/x//'g sample_\${SEGMENT}.fasta # To jest zbedne ?
        prep_own_vcf.py reference_\${SEGMENT}.fasta sample_\${SEGMENT}.fasta \$N \${vcf_input} \${vcf_output}
        bgzip \${vcf_output}; tabix \${vcf_output}.gz
        LISTA_VCF+=(\${vcf_output}.gz)
      done
      # na koniec concatujmey wszystkie vcf-y
     
      bcftools concat \${LISTA_VCF[@]}| bcftools sort --output-type z > ${sampleId}_final.vcf.gz
      tabix ${sampleId}_final.vcf.gz
   fi
   """
}
