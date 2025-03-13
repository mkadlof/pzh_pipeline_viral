#!/bin/bash


# Simplified script for running the BACTERIAL pipeline
# It is intended to be EXECUTED ON A100 machine
# The only parameters that do not have default values are:
# (I) directory with reads
# (II) sequencing platform
# all OTHER parameters ( localization od all the databases, images and modules) have PREDEFINED values that, if one chooses to, can be still modified (for testing purpose)

# localization of the main file with the pipeline and databases required to execute this pipeline
## Existing directories, for testing purpose can be change. but for production invariable
projectDir="/home/michall/git/pzh_pipeline_viral/"
external_databases_path="/mnt/raid/external_databases"
results_dir="./results" 

# docker images required to execute this pipeline
## Existing images, for testing purpose can be change, but for production invariable
main_image="pzh_pipeline_bacterial_main:latest"
prokka_image="staphb/prokka:latest"
alphafold_image="alphafold2:latest"

## Nextflow executor
profile="local"

# Parmaters related to resources available to the pipeline (max PER sample) if N samples are analyzed the pipeline will use at most N times more resuorces
# For testing purpose can be change, but for production invariable
threads=40

# Parameters without DEFAULTS that MUST be specified by a user
machine="" # Only Nanopore or Illumina
reads="" # Existing path


# Parameters with default values that depend on the sequencing platform
## For both platforms
genus="" # Salmonella Escherichia or Campylobacter can be empty
quality=""
min_number_of_reads="" 
min_median_quality=""
main_genus_value=""
kmerfinder_coverage=""
main_species_coverage=""
min_genome_length=""
unique_loci=""
contig_number=""
L50=""
final_coverage=""
min_coverage_ratio=""
min_coverage_value=""
## Nanopore-specific
model_medaka=""


# Usage function to display help
usage() {
    echo "Usage/Wywolanie: $0 --machine [Nanopore|Illumina] --reads PATH --projectDir PATH --external_databases_path PATH --main_image VALUE --prokka_image --alphafold_image VALUE[options]"
    echo "Required parameters/Parametry wymagane:"
    echo "  --machine VALUE                 Sequencing platform: Nanopore or Illumina"
    echo "                                   Platforma sekwencjonujaca uzyta do analizy. Mozliwe wartosci to Nanopore albo Illumina"
    echo "  --reads PATH                    Path to sequencing data with naming pattern for sequencing files. Use single quotes for this argument"
    echo "                                  Scieżka do katalogu z wynikami sekwencjonowania wraz z wzorcem nazewnictwa plików"
    echo "                                  Format plikow: fastq.gz, Przyklad: '/some/directory/*_R{1,2}.fastq.gz'"
    echo "  --projectDir PATH               Sciezka do katalogu z pobranym repozytorium"
    echo "                                  Directory with projects repository"
    echo "  --external_databases_path PATH  Sciezka do katalogu z pobranymi zewnetrznymi bazami"
    echo "                                  Directory with all databases used by the program"
    echo "  --main_image VALUE              Nazwa obrazu w formacie \"name:tag\" z obrazem zawierajacym programy uzywane przez pipeline"
    echo "                                  Name of the docker image with main program"
    echo "  --prokka_image VALUE            Nazwa obrazu w formacie \"name:tag\" z obrazem zawierajacym program prokka."
    echo "                                  Name of the docker image with prokka program"
    echo "  --alphafold_image VALUE         Nazwa obrazu w formacie \"name:tag\" z obrazem zawierajacym program alphafold"
    echo "                                  Name of the docker image with alphafold program"
    echo "Optional parameters:"
    echo "  --genus VALUE                   Genus of the bacteria that underwent sequencing"
    echo "                                  Nazwa rodzajowa bakterii podelgajacej sekwencjonowaiu"
    echo "                                  Akceptowane wartosci to Salmonella Escherichia lub Campylobacter"
    echo "                                  parametru nie trzeba podawac, program sam wykrywa rodzaj na podstawie odczytow"
    echo "  --threads VALUE                 Thread count (default: $threads)"
    echo "                                  Maksymalna ilosci CPU uzywana do analizy pojedycznej probki"
    echo "  --profile VALUE                 Nazwa profile zdefiniowanego w pliku konfiguaracyjnym nextflow z informacja o executor"
    echo "                                  Name of the profile specified in the nextflow configuration file."
    echo "                                  Available values: \"local\" and \"slurm\". Default value \"local\"."
    echo "  --all                           Display hidden parameters for advanced configuration"
    echo "                                  Wyswietl liste wszystkich parametrow modelu"
    echo "  -h, --help                      Show this help message"
}

# Full help
show_all_parameters() {
    echo "NOT YET IMPLEMENTED"
}


# Parse command-line options using GNU getopt
OPTS=$(getopt -o h --long projectDir:,profile:,external_databases_path:,results_dir:,main_image:,prokka_image:,alphafold_image:,threads:,machine:,reads:,genus:,quality:,min_number_of_reads:,min_median_quality:,main_genus_value:,kmerfinder_coverage:,main_species_coverage:,min_genome_length:,unique_loci:,contig_number:,L50:,final_coverage:,min_coverage_ratio:,min_coverage_value:,model_medaka:,all,help -- "$@")

eval set -- "$OPTS"

if [[ $# -eq 1 ]]; then
    echo "No parameters provided"
    usage
    exit 1
fi

while true; do
  case "$1" in
    --projectDir )
      projectDir="$2";
      shift 2
      ;;
    --profile)
      profile="$2"
      shift 2
      ;;
    --external_databases_path)
      external_databases_path="$2"
      shift 2
      ;;
    --results_dir)
      results_dir="$2"
      shift 2
      ;;
    --main_image )
      main_image="$2"; 
      shift 2 
      ;;
    --prokka_image )
      prokka_image="$2"; 
      shift 2 
      ;;
    --alphafold_image )
      alphafold_image="$2"; 
      shift 2 
      ;;
    --threads )
      threads="$2"; 
      shift 2 
      ;;
    --machine )
      machine="$2"; 
      shift 2 
      ;;
    --reads )
      reads="$2"; 
      shift 2 
      ;;
    --genus )
      genus="$2"; 
      shift 2 
      ;;
    --quality )
      quality="$2"; 
      shift 2 
      ;;
    --min_number_of_reads )
      min_number_of_reads="$2"; 
      shift 2 
      ;;
    --min_median_quality )
      min_median_quality="$2"; 
      shift 2 
      ;;
    --main_genus_value )
      main_genus_value="$2"; 
      shift 2 
      ;;
    --kmerfinder_coverage )
      kmerfinder_coverage="$2"; 
      shift 2 
      ;;
    --main_species_coverage )
      main_species_coverage="$2"; 
      shift 2 
      ;;
    --min_genome_length )
      min_genome_length="$2"; 
      shift 2 
      ;;
    --unique_loci )
      unique_loci="$2"; 
      shift 2 
      ;;
    --contig_number )
      contig_number="$2"; 
      shift 2 
      ;;
    --L50 )
      L50="$2"; 
      shift 2 
      ;;
    --final_coverage )
      final_coverage="$2"; 
      shift 2 
      ;;
    --min_coverage_ratio )
      min_coverage_ratio="$2"; 
      shift 2 
      ;;
    --min_coverage_value )
      min_coverage_value="$2"; 
      shift 2 
      ;;
    --model_medaka )
      model_medaka="$2"; 
      shift 2 
      ;;
    --all)
      show_all_parameters
      exit 0
      ;;
    -h|--help)
       usage
       exit 0
      ;;
    -- )
      shift; 
      break 
      ;;
    * )
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done



# Validate if all required parameters are set by a user
if [[ -z "$machine" || -z "$reads" ]]; then
    echo "Error: Missing required parameters."
    usage
    exit 1
fi

# Check if users provided one of the profiles understood by the pipeline
if [[ "$profile" != "slurm" && "$profile" != "local" ]]; then
    echo "User must provide one the available profiles: slurm or local"
    usage
    exit 1
fi

# Check if user provided correct sequencing platform and if so set default values for the main program
# Machine independent 
[[ -z "${min_coverage_ratio}" ]] && min_coverage_ratio=0.1
[[ -z "${min_coverage_value}" ]] && min_coverage_value=20
[[ -z "${main_genus_value}" ]] && main_genus_value=50
[[ -z "${kmerfinder_coverage}" ]] && kmerfinder_coverage=20
[[ -z "${main_species_coverage}" ]] && main_species_coverage=20
[[ -z "${min_genome_length}" ]] && min_genome_length=0.75
[[ -z "${L50}" ]] && L50=30000
[[ -z "${final_coverage}" ]] && final_coverage=20

if [[ "$machine" == "Illumina" ]]; then
        [[ -z "${quality}" ]] && quality=6
	[[ -z "${min_number_of_reads}" ]] && min_number_of_reads=50000
	[[ -z "${min_median_quality}" ]] && min_median_quality=10
	[[ -z "${unique_loci}" ]] && unique_loci=5
	[[ -z "${contig_number}" ]] && contig_number=1000
	[[ -z "${unique_loci}" ]] && unique_loci=5
elif [[ "$machine" == "Nanopore" ]]; then
        [[ -z "${quality}" ]] && quality=2
        [[ -z "${min_number_of_reads}" ]] && min_number_of_reads=10000
        [[ -z "${min_median_quality}" ]] && min_median_quality=5
        [[ -z "${contig_number}" ]] && contig_number=100
	[[ -z "${model_medaka}" ]] && model_medaka="r941_min_hac_g507"
	[[ -z "${unique_loci}" ]] && unique_loci=0
else
    echo "Error: Unsupported sequencinf platform: $machine. Supported values are: Nanopore, Illumina."
    exit 1
fi


# In a directory there must be at least one file that meet provided pattern
expanded_reads=$(eval ls ${reads} 2> /dev/null)

if [ $(echo "${expanded_reads}" | wc -w) -lt 1 ]; then
	echo "Error: No reads found in: ${reads}"
	exit 1
fi


# Tests for USER provided parameters 
# TO DO

echo "Running the bacterial pipeline..."
nextflow run ${projectDir}/nf_pipeline_bacterial.nf \
	     --results_dir ${results_dir} \
	     --genus ${genus} \
	     --reads "${reads}" \
	     --machine ${machine} \
	     --main_image ${main_image} \
	     --prokka_image ${prokka_image} \
	     --alphafold_image ${alphafold_image} \
	     --threads ${threads} \
	     --db_absolute_path_on_host ${external_databases_path} \
	     --min_coverage_ratio ${min_coverage_ratio} \
	     --min_coverage_value ${min_coverage_value} \
	     --quality ${quality} \
	     --min_number_of_reads ${min_number_of_reads} \
	     --min_median_quality ${min_median_quality} \
	     --main_genus_value ${main_genus_value} \
	     --kmerfinder_coverage ${kmerfinder_coverage} \
	     --main_species_coverage ${main_species_coverage} \
	     --min_genome_length ${min_genome_length} \
	     --unique_loci ${unique_loci} \
	     --contig_number ${contig_number} \
	     --L50 ${L50} \
	     --final_coverage ${final_coverage} \
	     --model_medaka ${model_medaka} \
	     -profile ${profile} \
	     -with-trace
