process varScan {
    tag "varScan:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    // publishDir "${params.results_dir}/${sampleId}/varscan", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(ref_genome)

    output:
    tuple val(sampleId), path('varscan.fa'), emit: fasta
    tuple val(sampleId), path('detected_variants_varscan.txt'), emit: pre_vcf

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch varscan.fa
      touch detected_variants_varscan.txt
    else
      # in case of SARS overlaping amplicons can have coverage twice params.max_depth so we adjust samtools accordingly
      MAX_DEPTH=`echo ${params.max_depth} | awk '{print int(\$0 * 1)}'` 
      samtools mpileup -B --max-depth \${MAX_DEPTH} \
                          --fasta-ref ${ref_genome} \
                          --min-BQ ${params.quality_snp} \
                          ${bam} >> ${bam}.mpileup

      varscan_qual=`echo "${params.quality_snp} - 1" | bc -l`
      java -jar /opt/varscan/VarScan.v2.4.6.jar pileup2cns ${bam}.mpileup \
                                                           --min-avg-qual \${varscan_qual} \
                                                           --p-value ${params.pval} \
                                                           --min-var-freq ${params.lower_ambig} \
                                                           --min-coverage ${params.min_cov} \
                                                           --variants \
                                                           --min-reads2 0 > detected_variants_varscan.txt

      parse_vcf_output_final.py detected_variants_varscan.txt ${params.upper_ambig} ${params.pval}

      bgzip --force detected_variants_varscan.vcf
      tabix detected_variants_varscan.vcf.gz

      qual=`echo ${params.pval} | awk '{print int(10*-log(\$1)/log(10))}'`

      bcftools norm --check-ref w \
                    --rm-dup all \
                    --fasta-ref ${ref_genome} \
                    detected_variants_varscan.vcf.gz | \
      bcftools norm --check-ref w \
                    --multiallelics -indels \
                    --fasta-ref ${ref_genome} | \
      bcftools filter \
                      --include "QUAL >= \${qual} && AF >= ${params.lower_ambig} && DP >= ${params.min_cov}" > detected_variants_varscan_final.vcf

      bgzip --force detected_variants_varscan_final.vcf
      tabix detected_variants_varscan_final.vcf.gz

      cat ${ref_genome} | bcftools consensus --mark-del X --samples - detected_variants_varscan_final.vcf.gz > varscan.fa
    fi
    """
}
