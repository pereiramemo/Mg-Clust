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
Usage: ./mg-clust_module-2.bash <options>
--help                          print this help
--assembly_file CHAR            input assembled metagenome (fasta format)
--bam_file CHAR                 bam input file (reads mapped to contigs)
--nslots NUM                    number of threads used (default 12)
--output_dir CHAR               directory to output generated data (default mg-clust_output-1)
--overwrite t|f                 overwrite previous folder if present (default f)
--sample_name CHAR              sample name used to name the files        
--train_file_name               train file name used to run FragGeneScan (default illumina_1)
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
  --assembly_file)
  if [[ -n "${2}" ]]; then
    ASSEMBLY_FILE="${2}"
    shift
  fi
  ;;
  --assembly_file=?*)
  ASSEMBLY_FILE="${1#*=}" 
  ;;
  --assembly_file=)
  printf 'Using default environment.\n' >&2
  ;;        
#############
  --bam_file)
  if [[ -n "${2}" ]]; then
    BAM_FILE="${2}"
    shift
  fi
  ;;
  --bam_file=?*)
  BAM_FILE="${1#*=}" 
  ;;
  --bam_file=)
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
#############
  --train_file_name)
  if [[ -n "${2}" ]]; then
    TRAIN_FILE_NAME="${2}"
    shift
  fi
  ;;
  --train_file_name=?*)
  TRAIN_FILE_NAME="${1#*=}" 
  ;;
  --train_file_name=)
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

if [[ ! -f "${ASSEMBLY_FILE}" ]]; then
  echo "input assembly file is not a real file"
  exit 1
fi  

if [[ ! -f "${BAM_FILE}" ]]; then
  echo "input bam file is not a real file"
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
    
    mkdir -p "${OUTPUT_DIR}"
    
    if [[ "$?" -ne "0" ]]; then
      echo "mkdir output directory ${OUTPUT_DIR} failed"
      exit 1
    fi
    
  elif [[  "${OVERWRITE}" == "f" ]]; then
    echo "${OUTPUT_DIR} already exists; use --overwrite t to overwrite"
    exit 0
  fi
  
elif [[ ! -d "${OUTPUT_DIR}" ]]; then

  mkdir -p "${OUTPUT_DIR}"
  
  if [[ "$?" -ne "0" ]]; then
    echo "mkdir output directory ${OUTPUT_DIR} failed"
    exit 1
  fi  
fi

###############################################################################
### 6. Predict ORFs
###############################################################################
 
"${fraggenescan}" \
-s "${ASSEMBLY_FILE}" \
-o "${OUTPUT_DIR}/${SAMPLE_NAME}_orfs" \
-w 0 \
--unordered \
-p "${NSLOTS}" \
-t "${TRAIN_FILE_NAME}" \
-r "${FRAGGENESCAN_TRAIN_FILE_DIR}" 

if [[ "$?" -ne 0 ]]; then
  echo "FragGeneScan failed"
  exit 1
fi  

###############################################################################
### 7. Create BED file
###############################################################################
  
FFN_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_orfs.ffn" 
BED_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_orfs.bed"

# create bed file
egrep ">" "${FFN_FILE}" | \
awk -v FS="_" -v OFS="\t" '{ 
  if ($0 ~ /^>/) {
  
    contig_prefix = $1
    sub(">", "", contig_prefix)
    contig_id = $2
    start = $3
    end = $4
    
    print contig_prefix"_"contig_id, start -1, end
  }
}' > "${BED_FILE}"

# Bed file is created zero-based for the start posistion 
# and one-based for the end position
# This is according to https://bedtools.readthedocs.io/en/latest/content/general-usage.html

if [[ "$?" -ne 0 ]]; then
  echo "awk command to create bed file failed"
  exit 1
fi  

###############################################################################
### 8. Get ORFs coverage
###############################################################################

# compute coverage: sorted markdup bam file
"${bedtools}" coverage \
-a "${BED_FILE}" \
-b "${BAM_FILE}" \
-mean > "${OUTPUT_DIR}/${SAMPLE_NAME}_orfs_tmp.cov"

if [[ "$?" -ne 0 ]]; then
  echo "bedtools to compute coverage failed"
  exit 1
fi  

###############################################################################
### 9. Fix ORFs start coordinates (back to one-based)
###############################################################################

awk 'BEGIN {OFS="\t"; FS="\t"} {

  contig_id = $1
  start_coord = $2 +1
  end_coord = $3
  abund = $4
  
  print contig_id, start_coord, end_coord, abund

}' "${OUTPUT_DIR}/${SAMPLE_NAME}_orfs_tmp.cov" > \
   "${OUTPUT_DIR}/${SAMPLE_NAME}_orfs.cov"

if [[ "$?" -ne 0 ]]; then
  echo "awk command to format start coords failed"
  exit 1
fi  

rm "${OUTPUT_DIR}/${SAMPLE_NAME}_orfs_tmp.cov"

if [[ "$?" -ne "0" ]]; then
  echo "removing ${SAMPLE_NAME}_orfs_tmp.cov failed"
  exit 1
fi  
  
###############################################################################
### 10. Exit
###############################################################################

echo "mg-clust_module-2.bash exited successfully"
exit 0
