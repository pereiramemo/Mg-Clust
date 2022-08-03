###############################################################################
### 1. Set env
###############################################################################

set -o pipefail 

source "/nfs/bin/mg-clust/mg-clust_conf.bash"

###############################################################################
### 2. Define help
###############################################################################

show_usage(){
  cat <<EOF
Usage: ./mg-clust_module-1.bash <options>
--help                          print this help
--assem_dir CHAR                directory with previously computed assemblies (format dirname/SAMPLE_NAME/SAMPLE_NAME.contigs.fa)
--assem_preset CHAR             MEGAHIT preset to generate assembly (default meta-sensitive)
--nslots NUM                    number of threads used (default 12)
--min_contig_length NUM         minimun length of contigs (smaller than this will be descarded; default 250) 
--output_dir CHAR               directory to output generated data (default mg-clust_output-1)
--overwrite t|f                 overwrite previous folder if present (default f)
--reads1 CHAR                   input R1 metagenome data (as fasta/q file)
--reads2 CHAR                   input R2 metagenome data (as fasta/q file)
--sample_name CHAR              sample name used to name the files        
EOF
}

###############################################################################
### 3. Parse input parameters
###############################################################################

while :; do
  case "${1}" in
    --help) # Call a "show_help" function to display a synopsis, then exit.
    show_usage
    exit 1;
    ;;
#############
  --assem_dir)
  if [[ -n "${2}" ]]; then
    ASSEM_DIR="${2}"
    shift
  fi
  ;;
  --assem_dir=?*)
  ASSEM_DIR="${1#*=}" 
  ;;
  --assem_dir=)
  printf 'Using default environment.\n' >&2
  ;;    
#############
  --assem_preset)
  if [[ -n "${2}" ]]; then
    ASSEM_PRESET="${2}"
    shift
  fi
  ;;
  --assem_preset=?*)
  ASSEM_PRESET="${1#*=}" 
  ;;
  --assem_preset=)
  printf 'Using default environment.\n' >&2
  ;;        
#############
  --reads1)
  if [[ -n "${2}" ]]; then
    R1="${2}"
    shift
  fi
  ;;
  --reads1=?*)
  R1="${1#*=}" 
  ;;
  --reads1=)
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --reads2)
  if [[ -n "${2}" ]]; then
    R2="${2}"
    shift
  fi
  ;;
  --reads2=?*)
  R2="${1#*=}" 
  ;;
  --reads2=)
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --nslots)
  if [[ -n "${2}" ]]; then
    NSLOTS="${2}"
    shift
  fi
  ;;
  --nslots=?*)
  NSLOTS="${1#*=}"
  ;;
  --nslots=)
  printf 'Using default environment.\n' >&2
  ;;
#############
  --min_contig_length)
  if [[ -n "${2}" ]]; then
    MIN_CONTIG_LENGTH="${2}"
    shift
  fi
  ;;
  --min_contig_length=?*)
  MIN_CONTIG_LENGTH="${1#*=}"
  ;;
  --min_contig_length=) 
  printf 'Using default environment.\n' >&2
  ;;
#############
  --output_dir)
  if [[ -n "${2}" ]]; then
    OUTPUT_DIR="${2}"
    shift
  fi
  ;;
  --output_dir=?*)
  OUTPUT_DIR="${1#*=}" 
  ;;
  --output_dir=)
  printf 'Using default environment.\n' >&2
  ;;
#############
  --overwrite)
  if [[ -n "${2}" ]]; then
    OVERWRITE="${2}"
    shift
  fi
  ;;
  --overwrite=?*)
  OVERWRITE="${1#*=}" 
  ;;
  --overwrite=)
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --sample_name)
  if [[ -n "${2}" ]]; then
    SAMPLE_NAME="${2}"
    shift
  fi
  ;;
  --sample_name=?*)
  SAMPLE_NAME="${1#*=}" 
  ;;
  --sample_name=)
  printf 'Using default environment.\n' >&2
  ;;  
############ End of all options.
  --)       
  shift
  break
  ;;
  -?*)
  printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
  ;;
  *) # Default case: If no more options, then break out of the loop.
  break
  esac
  shift
done  

###############################################################################
### 4. Check mandatory files
###############################################################################

if [[ ! -f "${R1}" ]]; then
  echo "read1 is not a real files"
  exit 1
fi  

if [[ ! -f "${R2}" ]]; then
  echo "read2 is not a real files"
  exit 1
fi  

###############################################################################
### 5. Create output folder
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then
  if [[ "${OVERWRITE}" == "t" ]]; then
    
    rm -r "${OUTPUT_DIR}"
    
    if [[ $? -ne "0" ]]; then
      echo "rm -r output directory ${OUTPUT_DIR} failed"
      exit 1
    fi
    
  elif [[  "${OVERWRITE}" == "f" ]]; then
    echo "${OUTPUT_DIR} already exists; use --overwrite t to overwrite"
    exit 0
  fi
  
fi

###############################################################################
### 6. De novo assembly
###############################################################################

if [[ -z "${ASSEM_DIR}" ]]; then

  # ${OUTPUT_DIR} is created by megahit
  "${megahit}" \
  --num-cpu-threads "${NSLOTS}" \
  -1 "${R1}" \
  -2 "${R2}" \
  --presets ${ASSEM_PRESET} \
  --min-contig-len "${MIN_CONTIG_LENGTH}" \
  --out-prefix "${SAMPLE_NAME}" \
  --out-dir "${OUTPUT_DIR}" 

  if [[ "$?" -ne 0 ]]; then
    echo "megahit failed"
    exit 1
  fi  

  ASSEMBLY_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.contigs.fa"
  
else 

  mkdir "${OUTPUT_DIR}"
   if [[ "$?" -ne "0" ]]; then
    echo "mkdir ${OUTPUT_DIR} failed"
    exit 1
  fi 
  
  ASSEMBLY_FILE="${ASSEM_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.contigs.fa"

fi

###############################################################################
### 7. Map short reads
###############################################################################

ASSEMBLY_FILE_NSEQ=$(egrep -c ">" "${ASSEMBLY_FILE}")

if [[ "${ASSEMBLY_FILE_NSEQ}" -lt 5 ]]; then
  echo "Not enough assembled sequences to continue: ${ASSEMBLY_FILE_NSEQ}"
  exit 0
fi  

# index
"${bwa}" index "${ASSEMBLY_FILE}"
  
if [[ $? -ne 0 ]]; then
  echo "bwa index failed"
  exit 1
fi  
 
# map
"${bwa}" mem -M -t "${NSLOTS}" "${ASSEMBLY_FILE}" "${R1}" "${R2}" > \
"${OUTPUT_DIR}/${SAMPLE_NAME}.sam"
  
if [[ $? -ne 0 ]]; then
  echo "bwa mem failed"
  exit 1
fi  
  
# convert to bam and filter high-quality primary alignments
# flags taken from https://broadinstitute.github.io/picard/explain-flags.html
# secondary alignments are moved; however this makes very little difference
# ORFs coverage using F4 vs F260 flags had a MSE=0.00486  and Pearson cor 0.997 (same toydataset)

"${samtools}" view \
-@ "${NSLOTS}" \
-q 10 \
-F 260 \
-b "${OUTPUT_DIR}/${SAMPLE_NAME}.sam" > "${OUTPUT_DIR}/${SAMPLE_NAME}.bam"
  
if [[ $? -ne 0 ]]; then
  echo "samtools convert to bam failed"
  exit 1
fi  
  
# sort 
"${samtools}" \
sort -@ "${NSLOTS}" \
"${OUTPUT_DIR}/${SAMPLE_NAME}.bam" > "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"
  
if [[ $? -ne 0 ]]; then
  echo "samtools sort failed"
  exit 1
fi
  
# index merged data
"${samtools}" index -@ "${NSLOTS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"
  
if [[ $? -ne 0 ]]; then
  echo "samtools inex failed"
  exit 1
fi

# remove duplicates (defined as originating from a single fragment of DNA)
mkdir "${OUTPUT_DIR}/tmp"

if [[ $? -ne 0 ]]; then
  echo "mkdir ${OUTPUT_DIR}/tmp failed" 
  exit 1
fi  

java -jar "${picard}" MarkDuplicates \
INPUT="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.bam" \
METRICS_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.metrics.txt" \
REMOVE_DUPLICATES=TRUE \
ASSUME_SORTED=TRUE \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 \
TMP_DIR="${OUTPUT_DIR}/tmp"
  
if [[ $? -ne 0 ]]; then
  echo "picard failed"
  exit 1
fi
  
###############################################################################
### 8. Clean
###############################################################################
  
rm -r "${OUTPUT_DIR}/${SAMPLE_NAME}.sam" \
      "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" \
      "${OUTPUT_DIR}/tmp"

if [[ $? -ne "0" ]]; then
  echo "removing intermediate mapping files failed"
  exit 1
fi  
      
if [[ -z "${ASSEM_DIR}" ]]; then
  rm -r "${OUTPUT_DIR}/intermediate_contigs"
    
  if [[ $? -ne 0 ]]; then
    echo "removing assembly intermediate files failed"
    exit 1
  fi
  
fi
  
###############################################################################
### 9. Exit
###############################################################################

echo "mg-clust_module-1.bash exited successfully"
exit 0
