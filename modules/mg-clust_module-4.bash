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
Usage: ./mg-clust_module-4.bash <options>
--help                          print this help
--clust_thres NUM               clustering threshold (default 0.7)
--cov_table CHAR                ORFs' coverage table
--min_opu_occup NUM             minimum OPU occupancy (smaller than this will be discarded; default 2)
--nslots NUM                    number of threads used (default 12)
--orfs_db CHAR                  mmseqs orfs db
--output_dir CHAR               directory to output generated data (default mg-clust_output-2)
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
  --cov_table)
  if [[ -n "${2}" ]]; then
    COV_TABLE="${2}"
    shift
  fi
  ;;
  --cov_table=?*)
  COV_TABLE="${1#*=}" 
  ;;
  --cov_table=)
  printf 'Using default environment.\n' >&2
  ;;     
#############
  --clust_thres)
  if [[ -n "${2}" ]]; then
    CLUST_THRES="${2}"
    shift
  fi
  ;;
  --clust_thres=?*)
  CLUST_THRES="${1#*=}" 
  ;;
  --clust_thres=)
  printf 'Using default environment.\n' >&2
  ;;
#############
  --min_opu_occup)
  if [[ -n "${2}" ]]; then
    MIN_OPU_OCCUP="${2}"
    shift
  fi
  ;;
  --min_opu_occup=?*)
  MIN_OPU_OCCUP="${1#*=}" 
  ;;
  --min_opu_occup=)
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
  --orfs_db)
  if [[ -n "${2}" ]]; then
    ORFS_DB="${2}"
    shift
  fi
  ;;
  --orfs_db=?*)
  ORFS_DB="${1#*=}" 
  ;;
  --orfs_db=)
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
### 4. Run orfs clustering
###############################################################################

CLUST_DIR="${OUTPUT_DIR}/clust_orfs_id${CLUST_THRES/0.}"

mkdir "${CLUST_DIR}"
if [[ "$?" -ne "0" ]]; then 
  echo "mkdir ${CLUST_DIR} failed"
  exit 1
fi  

if [[ -d "${CLUST_DIR}/tmp" ]]; then
  rm -r "${CLUST_DIR}/tmp"
  if [[ "$?" -ne "0" ]]; then 
    echo "rm -r ${CLUST_DIR} failed"
    exit 1
  fi  
fi    
    
"${mmseqs}" cluster \
"${ORFS_DB}" \
"${CLUST_DIR}/orfs_clust_id${CLUST_THRES/0.}" \
"${CLUST_DIR}/tmp" \
--min-seq-id "${CLUST_THRES}" \
--threads "${NSLOTS}" \
--cov-mode 0 \
-c 0.85

if [[ "$?" -ne "0" ]]; then
  echo "mmseqs2 cluster failed"
  exit 1
fi  

rm -r "${CLUST_DIR}/tmp"

if [[ "$?" -ne "0" ]]; then
  echo "rm ${CLUST_DIR}/tmp failed"
  exit 1
fi  

###############################################################################
### 5. Convert to tab tables
###############################################################################

CLUST_TABLE="${CLUST_DIR}/orfs_clust_id${CLUST_THRES/0.}.tsv"

"${mmseqs}" createtsv \
"${ORFS_DB}" \
"${ORFS_DB}" \
"${CLUST_DIR}/orfs_clust_id${CLUST_THRES/0.}" \
"${CLUST_TABLE}"

if [[ "$?" -ne "0" ]]; then
  echo "mmseqs createtsv failed"
  exit 1
fi  

###############################################################################
### 6. Cross tables
###############################################################################

CLUST2ABUND_TABLE="${CLUST_DIR}/orfs_clust2abund_id${CLUST_THRES/0.}.tsv"

awk -v FS="\t" -v OFS="\t" -v outdir="${CLUST_DIR}" -v id="${CLUST_THRES/0.}" \
'{
  
  if (NR == FNR) {
  
    clust_id=gensub("_[+,-]$","","g",$1) 
    seq_id=gensub("_[+,-]$","","g",$2)
    
    if ( array_clust[seq_id] == "") {
   
      array_clust[seq_id]=clust_id
      
    } else {
    
      print $1 > outdir"/duplicated_"id".list" 
      
      # seq_ids should not be duplicated unless there is an ORF in the exact 
      # same location in the opposite strand
      # In such a case, we only analyze the first ORF that appears in the table
      
    }
    
    next
  }   
      
  if (array_clust[$1] != "") { 
    
    clust_id=array_clust[$1]
    seq_id = $1
    abund = $2
    sample=gensub("(.*)-k[0-9]+_[0-9]+.*","\\1","g",seq_id)

    print sample, clust_id, seq_id, abund
    
  } else {
  
    print $1 > outdir"/not_found_"id".list"
    
    # seq_ids filtered by length (module 3) will not be found
    
  }  
     
}' "${CLUST_TABLE}" "${COV_TABLE}" > "${CLUST2ABUND_TABLE}"

if [[ "$?" -ne "0" ]]; then
  echo "map cluster to abund failed (awk command)"
  exit 1
fi  

###############################################################################
### 7. Create workables: collapsed and filtered
###############################################################################

PARENT_DIR=$(dirname "${OUTPUT_DIR}")
CLUST2ABUND_TABLE_WORKABLE="${PARENT_DIR}/orfs_clust2abund_id${CLUST_THRES/0.}_workable.tsv"

# We use gawk to facilitate multidimensional arrays
gawk -v FS="\t" -v OFS="\t" -v min_opu_occup="${MIN_OPU_OCCUP}" \
'{
  
  sample=$1
  opu_id=$2
  abund=$4
  
  #print sample, opu_id, abund
  
  clust2abund_array[sample][opu_id] = clust2abund_array[sample][opu_id] + abund  
  clust2occup_array[opu_id][sample] =  1
    
} END {

  for (sample in clust2abund_array) { 
    for (opu_id in clust2abund_array[sample]) {
    
      occup_total = length(clust2occup_array[opu_id])
      abund_total = clust2abund_array[sample][opu_id]
    
      if (occup_total >= min_opu_occup) {
    
        print sample,opu_id,abund_total
      }
    }    
  }

}' "${CLUST2ABUND_TABLE}" > "${CLUST2ABUND_TABLE_WORKABLE}"

if [[ "$?" -ne "0" ]]; then
  echo "collapse and filter opus failed (awk command)"
  exit 1
fi  

###############################################################################
### 8. Exit
###############################################################################

echo "mg-clust_module-4.bash exited successfully"
exit 0
