process nextclade {
    tag "nextclade:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "nextstrain_lineage.csv"
    container  = params.main_image
    containerOptions "--volume ${params.external_databases_path}:/home/external_databases/"
    cpus 1
    memory "20 GB"
    input:
    tuple val(sampleId), path('output_consensus_masked_SV.fa'), path(ref_genome_with_index), val(QC_status), val(SAMPLE_SUBTYPE)

    output:
    tuple val(sampleId), path("nextstrain_lineage.json"), emit: json
    tuple val(sampleId), path('dummy.fasta'), val(QC_status), emit:to_modeller

    script:
    """
    set -x
    touch dummy.fasta
    
 # Default json file
    generate_empty_json() {
        local segment="\$1"
        cat <<EOF > "nextstrain_lineage_\${segment}.json"
{
    "status": "nie",
    "database_name": "Nextclade",
    "error_message": "Nextclade was skip for ${SAMPLE_SUBTYPE} in sample ${sampleId} for segment \${segment}",
    "sequence_source": "\${segment}"
}
EOF
    }

    # we split final genome into segments
    if [ ${QC_status} == "nie" ]; then
        generate_empty_json "HA"
        generate_empty_json "NA"
        jq -s "." nextstrain_lineage_HA.json nextstrain_lineage_NA.json > nextstrain_lineage.json
    else
        cat output_consensus_masked_SV.fa | awk '{if (substr(\$0, 1, 1)==">") { new_name=\$0; gsub("\\\\.", "_", new_name); gsub("/", "_", new_name);  gsub("_SV", "", new_name);  filename=("sample_"substr(new_name,2) ".fasta"); print \$0 >> filename } else {print toupper(\$0)  >> filename}}'
        KNOWN='H1N1 H3N2 Yamagata Victoria'
        if [[ \${KNOWN[@]} =~ ${SAMPLE_SUBTYPE} ]]; then
            nextclade run \
                --input-dataset /home/external_databases/nextclade/${SAMPLE_SUBTYPE}_HA.zip \
                --output-csv nextclade_lineage_HA.csv \
                --output-all nextclade_linages_HA \
                sample_chr4_HA.fasta
            if [[ -s nextstrain_lineage_HA.csv ]]; then
                parse_nextclade_output_csv2json.py nextstrain_lineage_HA.csv nextclade_lineages_HA/nextclade.auspice.json nextstrain_lineage_HA.json HA
            fi
            # for all this subtypes, save, Yamagata, we can also analyze NA
            if [ ${SAMPLE_SUBTYPE} != 'Yamagata' ]; then
                nextclade run --input-dataset /home/external_databases/nextclade/${SAMPLE_SUBTYPE}_NA.zip \
                    --output-csv nextclade_lineage_NA.csv \
                    --output-all nextclade_linages_NA \
                    sample_chr6_NA.fasta
                if [[ -s nextstrain_lineage_NA.csv ]]; then
                    parse_nextclade_output_csv2json.py nextstrain_lineage_NA.csv nextclade_lineages_NA/nextclade.auspice.json nextstrain_lineage_NA.json NA
                fi
            fi
        fi

        # If files where not created by Nextclade then create empty ones.
        if [[ ! -f "nextstrain_lineage_HA.json" ]]; then
            generate_empty_json "HA"
        fi

        if [[ ! -f "nextstrain_lineage_NA.json" ]]; then
            generate_empty_json "NA"
        fi

        # Combine two jsons in one
        jq -s "." nextstrain_lineage_HA.json nextstrain_lineage_NA.json > nextstrain_lineage.json
    fi
    """
}
