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
Usage: ./mg-clust_module-3.bash <options>
--help                          print this help
--input_dir CHAR                input dir where to find the ORF fasta files and coverage tables
--nslots NUM                    number of threads used (default 12)
--min_orf_length NUM            minimum length of ORFs (amino acids); ORFs shorter than this will be discarded (default 60)
--output_dir CHAR               directory to output generated data
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
  --input_dir)
  if [[ -n "${2}" ]]; then
    INPUT_DIR="${2}"
    shift
  fi
  ;;
  --input_dir=?*)
  INPUT_DIR="${1#*=}" 
  ;;
  --input_dir=)
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
  --min_orf_length)
  if [[ -n "${2}" ]]; then
    MIN_ORF_LENGTH="${2}"
    shift
  fi
  ;;
  --min_orf_length=?*)
  MIN_ORF_LENGTH="${1#*=}" 
  ;;
  --orf_min_length=)
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
### 4. Concat ORF files
###############################################################################

CONCAT_ORFS="${OUTPUT_DIR}/orfs.faa"

ls "${INPUT_DIR}/"*"/"*"_orfs.faa" | \
while read LINE; do 

  SAMPLE=$(basename "${LINE}" _orfs.faa);
  cat "${LINE}" | sed "s/^>/>${SAMPLE}-/" 
  
  if [[ "$?" -ne "0" ]]; then
    echo "concat orf file ${LINE} failed"
    exit 1
  fi
  
done > "${CONCAT_ORFS}"

###############################################################################
### 5. Filter ORFs by length
###############################################################################

CONCAT_ORFS_FILTERED=${CONCAT_ORFS/.faa/_filtered.faa}

if [[ ! -f "${ORFS_FILTERED}" ]]; then

  "${bbduk}" \
  in="${CONCAT_ORFS}" \
  out="${CONCAT_ORFS_FILTERED}" \
  overwrite=t \
  threads="${NSLOTS}" \
  minlength="${MIN_ORF_LENGTH}" \
  amino=t

fi  
  
if [[ "$?" -ne "0" ]]; then
  echo "bbduk failed"
  exit 1
fi
  
###############################################################################
### 6. Create mmseqs database
###############################################################################

ORFS_DB=${CONCAT_ORFS_FILTERED/.faa/_db}

if [[ ! -f "${ORFS_DB}" ]]; then
  
  "${mmseqs}" createdb \
  "${CONCAT_ORFS_FILTERED}" \
  "${ORFS_DB}" \
  --dbtype 1 
  
  if [[ "$?" -ne "0" ]]; then
    echo "mmseqs createdb failed"
    exit 1
  fi

fi

###############################################################################
### 7. Concat and format coverage table
###############################################################################

COV_TABLE="${OUTPUT_DIR}/orfs_cov.tsv"

ls "${INPUT_DIR}/"*"/"*"_orfs.cov" |
while read LINE; do

  SAMPLE=$(basename "${LINE}" _orfs.cov)
  awk -v S="${SAMPLE}" -v OFS="\t" -v FS="\t" '{
  
    print S"-"$1"_"$2"_"$3,$4
  
  }' "${LINE}" 
  
  if [[ $? -ne "0" ]]; then
    echo "concat and format coverage table failed (awk command)"
    exit 1
  fi  

done > "${COV_TABLE}"

###############################################################################
### 8. Concat and format coverage table
###############################################################################

echo "mg-clust_module-3.bash exited successfully"
exit 0
