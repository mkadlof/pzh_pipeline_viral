process medaka_varscan_integration_first_round {
    // This module integrates medaka vcf and SPECIFIC parts of varscan prediction
    // i.e. SNPs with over 80% usage (according to varscan) that are not predicted as valid SNPs according to medaka
    // due to complex underlying genotype of a region
    tag "programs integration:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz'), path('medaka_annotated_filtered.vcf.gz.tbi'), path('genome.fasta'), val(QC_status),  path('detected_variants_varscan.txt')

    output:
    tuple val(sampleId), path('medaka_and_varscan_final.vcf.gz'), path('medaka_and_varscan_final.vcf.gz.tbi'), val(QC_status), emit: vcf
    tuple val(sampleId),  path('genome.fasta'),  val(QC_status), emit: reference_genome

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch medaka_and_varscan_final.vcf.gz
      touch medaka_and_varscan_final.vcf.gz.tbi
    else
        MINIMUM_USAGE=80
        awk -v usage="\${MINIMUM_USAGE}" '{if (int(\$7) > usage) print \$0; else if (\$1 ~ /Chrom/) print \$0}' detected_variants_varscan.txt >> detected_variants_varscan_part1_filtered.txt
        
        merge_varscan_with_medaka_final_INFL_first_round.py medaka_annotated_filtered.vcf.gz detected_variants_varscan_part1_filtered.txt  medaka_and_varscan.vcf
  
        # QUAL AND DP are insignificant becuase we filter mutations during medaka step
        bcftools sort medaka_and_varscan.vcf |  bcftools norm -c w -d all -f genome.fasta | bcftools norm -c w -m -indels -f genome.fasta | bcftools filter -O z -o medaka_and_varscan_final.vcf.gz -i "QUAL >= 0 && INFO/DP >= 1"

        tabix medaka_and_varscan_final.vcf.gz
    


    fi
    """
}


process medaka_varscan_integration_second_round {
    // This module integrates medaka vcf and SPECIFIC parts of varscan prediction
    // In second round varscan is used to introduce ambigous nucleotides
    tag "programs integration:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz'), path('medaka_annotated_filtered.vcf.gz.tbi'), path('medaka_annotated.vcf.gz'), path('medaka_annotated.vcf.gz.tbi'), path('genome.fasta'), val(QC_status),  path('detected_variants_varscan.txt')

    output:
    tuple val(sampleId), path('genome.fasta'), val(QC_status), path('medaka_and_varscan_final.vcf.gz'), path('medaka_and_varscan_final.vcf.gz.tbi'), emit: vcf
    tuple val(sampleId),  path('genome.fasta'), emit: reference_genome

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch medaka_and_varscan_final.vcf.gz
      touch medaka_and_varscan_final.vcf.gz.tbi
    else
      # This script requires two vcf from medaka filtered and unfiltered
     
      merge_varscan_with_medaka_final_INFL.py medaka_annotated_filtered.vcf.gz detected_variants_varscan.txt ${params.upper_ambig} medaka_and_varscan.vcf medaka_annotated.vcf.gz ${params.min_cov}

      
      # From final VCF write only non-frameshift mutations and frameshift mutations that DP is twice as much as the masking threshold
      # This was origianlly done only for SARS but we introduce it for all organisms
      MIN_COV_FRAMESHIFT=`echo "${params.min_cov} * 2" | bc -l `
      cat medaka_and_varscan.vcf | awk -v COV="\${MIN_COV_FRAMESHIFT}"  '{split(\$8, a, ";"); split(a[1], DP, "="); if (substr(\$0, 1, 1)=="#") {print \$0} else {if ((length(\$5) - length(\$4)) % 3 == 0)  {print \$0} else if (( (length(\$5) - length(\$4) ) % 3 != 0 && DP[2] > COV) )  {print \$0} }}' > non_frameshift.vcf

      # QUAL AND DP are insignificant becuase we filter mutations during medaka step
      bcftools sort non_frameshift.vcf |  bcftools norm -c w -d all -f genome.fasta | bcftools norm -c w -m -indels -f genome.fasta | bcftools filter -O z -o medaka_and_varscan_final.vcf.gz -i "QUAL >= 0 && INFO/DP >= 1"

      tabix medaka_and_varscan_final.vcf.gz 
    fi
    """
}
