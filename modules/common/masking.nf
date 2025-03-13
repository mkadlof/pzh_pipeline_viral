process masking {
    tag "masking:${sampleId}"
    container  = params.main_image
    input:
    tuple val(sampleId), path(bam), path(bai), path(primers), path(pairs), val(QC_status)
    cpus 1
    memory "20 GB"
    output:
    tuple val(sampleId), path('ivar_trimmed_all.bam'), val(QC_status)

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch ivar_trimmed_all.bam
    else
      length=`echo "${params.length} - 40" | bc -l`
      ivar trim -i ${bam} \
               -b ${primers} \
               -m \${length} \
               -f ${pairs} \
               -q ${params.quality_initial} \
               -e \
               -p ivar_trimmed_all 
    fi
    """
}

process masking_nanopore {
    // Ten modult przyjmuje poza glownym kanalem 
    // dwa dodatkowe. Pierwszy ktory zostanie zmapowany na cmienna tolerance jest parametrem samtools ampliconclip
    // drugi zostanie zmapowany na wartosc skip_trimmin aby dalo sie w przyszlosci
    // omisjac krok trimmingu gdy mamy do czynienia z sekwencjonowaniem metagenomicznym wirusow jak w EQA2024
    tag "masking:${sampleId}"
    container  = params.main_image
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path(bam), path(bai), path(genome), path(primers), val(QC_status)
    val(tolerance)
    val(skip_trimming)

    output:
    tuple val(sampleId), path('trimmed.bam'), path('trimmed.bam.bai'),  val(QC_status), path(genome), emit: bam_and_genome
    tuple val(sampleId), path('trimmed.bam'), path('trimmed.bam.bai'),  val(QC_status), emit: bam_only

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch trimmed.bam
      touch trimmed.bam.bai
    else
      if [ ${skip_trimming} -eq 1 ]; then
        samtools sort -o trimmed.bam $bam
        samtools index trimmed.bam
      else
        samtools ampliconclip --filter-len 1 --both-ends -b ${primers} --tolerance ${tolerance} -o forvariants_initial.bam -O bam ${bam} 
        samtools sort -o trimmed.bam forvariants_initial.bam
        samtools index trimmed.bam

      fi
    fi
    """
}

