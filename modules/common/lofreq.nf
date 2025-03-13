process lofreq {
    tag "lofreq:${sampleId}"
    container  = params.main_image
    // publishDir "${params.results_dir}/${sampleId}/lofreq", mode: 'copy'
    cpus { params.threads > 10 ? 10 : params.threads }
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(ref_genome)

    output:
    tuple val(sampleId), path('lofreq.fa')

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch lofreq.fa
    else
      # this module recieves only genome in fasta format 
      samtools faidx $ref_genome

      lofreq call-parallel --pp-threads ${task.cpus} \
                           --ref ${ref_genome} \
                           --max-depth ${params.max_depth} \
                           --min-cov ${params.min_cov} \
                           --call-indels \
                           --out detected_variants_lofreq.vcf \
                           ${bam}
    
      cat detected_variants_lofreq.vcf | \
          bcftools norm --check-ref w \
                        --rm-dup all \
                        --fasta-ref ${ref_genome} | \
          bcftools norm --check-ref w \
                         --multiallelics -indels \
                         --fasta-ref ${ref_genome} > detected_variants_lofreq_fix.vcf

      qual=`echo ${params.pval} | awk '{print int(10*-log(\$1)/log(10))}'`
      cat detected_variants_lofreq_fix.vcf | \
          awk -v qual=\${qual} \
              -v min_cov=${params.min_cov} \
              -v upper_ambig=${params.upper_ambig} \
              -v lower_ambig=${params.lower_ambig} \
              '{if(substr(\$1, 1, 1) == "#") print \$0; else split(\$8,a, ";"); split(a[4], b, "="); split(b[2], c, ","); split(a[1], DP, "="); split(a[2], AF, "="); if( ((c[3] + c[4]) / (c[3] + c[4] + c[1] + c[2])  > upper_ambig &&  \$6 >= qual  && DP[2] >= min_cov) || (a[5] == "INDEL" && AF[2] >= lower_ambig && \$6 >= qual  && DP[2] >= min_cov)) print \$0}' > detected_variants_lofreq_fix_high.vcf
    
      bgzip --force detected_variants_lofreq_fix_high.vcf
      tabix detected_variants_lofreq_fix_high.vcf.gz
    
      cat detected_variants_lofreq_fix.vcf | \
          awk -v qual=\${qual} \
              -v min_cov=${params.min_cov} \
              -v upper_ambig=${params.upper_ambig} \
              -v lower_ambig=${params.lower_ambig} \
              '{if(substr(\$1, 1, 1) == "#") print \$0; else split(\$8,a, ";"); split(a[4], b, "="); split(b[2], c, ","); split(a[1], DP, "="); if( (c[3] + c[4]) / (c[3] + c[4] + c[1] + c[2])  <= upper_ambig &&  \$6 >= qual  && DP[2] >= min_cov && (c[3] + c[4]) / (c[3] + c[4] + c[1] + c[2])  >= lower_ambig && (a[5] != "INDEL")) print \$0}' > tmp_lofreq.vcf
    
      introduce_amb_2_vcf.py tmp_lofreq.vcf \
                             detected_variants_lofreq_fix_ambig.vcf
   
      bgzip --force detected_variants_lofreq_fix_ambig.vcf
      tabix detected_variants_lofreq_fix_ambig.vcf.gz
    
      bcftools concat detected_variants_lofreq_fix_high.vcf.gz \
                      detected_variants_lofreq_fix_ambig.vcf.gz | \
                      bcftools sort --output-type z > detected_variants_lofreq_final.vcf.gz
    
      tabix detected_variants_lofreq_final.vcf.gz
    
      cat ${ref_genome} | \
          bcftools consensus --mark-del X --samples - \
          detected_variants_lofreq_final.vcf.gz > lofreq.fa
    fi # koniec if-a na zle QC 
    """
}





