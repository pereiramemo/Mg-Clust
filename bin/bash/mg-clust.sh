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
Usage: ./mg-clust.bash <options>
--help                          print this help
--assem_dir CHAR                directory with previously computed assemblies (format dirname/SAMPLE_NAME/SAMPLE_NAME.contigs.fa)
--assem_preset CHAR             MEGAHIT preset to generate assembly (default meta-sensitive)
--compress t|f                  compress all output data (default f)
--clean t|f                     clean up intermediate data (default f)
--input_dir CHAR                directory of input metagenomes
--logs_file CHAR                file name to save parallel logs
--nslots NUM                    number of threads used (default 12)
--njobs NUM                     number of jobs to run in parallel (each job with nslots) (default 3)
--min_contig_length NUM         minimum length of contigs (smaller than this will be discarded; default 250) 
--min_opu_occup NUM             minimum OPU occupancy (smaller than this will be discarded; default 2)
--min_orf_length NUM            minimum length of ORFs (amino acids); ORFs shorter than this will be discarded (default 60)
--output_dir CHAR               directory to output generated data (default mg-clust_output)
--overwrite t|f                 overwrite previous folder if present (default f)
--reads1_suffix CHAR            suffix of R1 reads
--reads2_suffix CHAR            suffix of R2 reads
--run_module_1 t|f              run the first processing module (assemble and map reads; this module will fail if folder output-1 exists; default t)
--run_module_2 t|f              run the second processing module (predict ORFs and compute ORFs coverage; this module will fail if folder output-2 exists; default t)
--run_module_3 t|f              run the third processing module (concatenate data and create ORFs db; this module will fail if folder output-3 exists; default t)
--run_module_4 t|f              run the fourth processing module (cluster ORFs and compute clusters abundance; folder output-4 will be kept if present; default t)
--servers CHAR,CHAR             comma separated list of servers to run mg-clust
--train_file_name               train file name used to run FragGeneScanRs (default illumina_1)
--thres_range NUM,NUM           minimum and maximum clustering thresholds separated by comma (default 0.7,0.9)
--thres_step NUM                threshold sequence step (default 0.1)
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
  --compress)
  if [[ -n "${2}" ]]; then
    COMPRESS="${2}"
    shift
  fi
  ;;
  --compress=?*)
  COMPRESS="${1#*=}" 
  ;;
  --compress=)
  printf 'Using default environment.\n' >&2
  ;;    
#############
  --clean)
  if [[ -n "${2}" ]]; then
    CLEAN="${2}"
    shift
  fi
  ;;
  --clean=?*)
  CLEAN="${1#*=}" 
  ;;
  --clean=)
  printf 'Using default environment.\n' >&2
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
  --logs_file)
  if [[ -n "${2}" ]]; then
    LOGS_FILE="${2}"
    shift
  fi
  ;;
  --logs_file=?*)
  LOGS_FILE="${1#*=}"
  ;;
  --logs_file=) 
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
  --min_orf_length)
  if [[ -n "${2}" ]]; then
    MIN_ORF_LENGTH="${2}"
    shift
  fi
  ;;
  --min_orf_length=?*)
  MIN_ORF_LENGTH="${1#*=}" 
  ;;
  --min_orf_length=)
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
  --njobs)
  if [[ -n "${2}" ]]; then
    NJOBS="${2}"
    shift
  fi
  ;;
  --njobs=?*)
  NJOBS="${1#*=}"
  ;;
  --njobs=)
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
  --reads1_suffix)
  if [[ -n "${2}" ]]; then
    R1_SUFFIX="${2}"
    shift
  fi
  ;;
  --reads1_suffix=?*)
  R1_SUFFIX="${1#*=}" 
  ;;
  --reads1_suffix=)
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --reads2_suffix)
  if [[ -n "${2}" ]]; then
    R2_SUFFIX="${2}"
    shift
  fi
  ;;
  --reads2_suffix=?*)
  R2_SUFFIX="${1#*=}" 
  ;;
  --reads2=)
  printf 'Using default environment.\n' >&2
  ;;
#############
  --run_module_1)
  if [[ -n "${2}" ]]; then
    RUN_MODULE_1="${2}"
    shift
  fi
  ;;
  --run_module_1=?*)
  RUN_MODULE_1="${1#*=}" 
  ;;
  --run_module_1=)
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --run_module_2)
  if [[ -n "${2}" ]]; then
    RUN_MODULE_2="${2}"
    shift
  fi
  ;;
  --run_module_2=?*)
  RUN_MODULE_2="${1#*=}" 
  ;;
  --run_module_2=)
  printf 'Using default environment.\n' >&2
  ;;    
#############
  --run_module_3)
  if [[ -n "${2}" ]]; then
    RUN_MODULE_3="${2}"
    shift
  fi
  ;;
  --run_module_3=?*)
  RUN_MODULE_3="${1#*=}" 
  ;;
  --run_module_3=)
  printf 'Using default environment.\n' >&2
  ;;   
#############
  --run_module_4)
  if [[ -n "${2}" ]]; then
    RUN_MODULE_4="${2}"
    shift
  fi
  ;;
  --run_module_4=?*)
  RUN_MODULE_4="${1#*=}" 
  ;;
  --run_module_4=)
  printf 'Using default environment.\n' >&2
  ;;    
#############
  --servers)
  if [[ -n "${2}" ]]; then
    SERVERS="${2}"
    shift
  fi
  ;;
  --servers=?*)
  SERVERS="${1#*=}" 
  ;;
  --servers=)
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
#############
  --thres_range)
  if [[ -n "${2}" ]]; then
    THRES_RANGE="${2}"
    shift
  fi
  ;;
  --thres_range=?*)
  THRES_RANGE="${1#*=}" 
  ;;
  --thres_range=)
  printf 'Using default environment.\n' >&2
  ;;    
#############
  --thres_step)
  if [[ -n "${2}" ]]; then
    THRES_STEP="${2}"
    shift
  fi
  ;;
  --thres_step=?*)
  THRES_STEP="${1#*=}" 
  ;;
  --thres_step=)
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
### 4. Define defaults
###############################################################################

if [[ -z "${ASSEM_PRESET}" ]]; then
  ASSEM_PRESET="meta-sensitive"
fi

if [[ -z "${COMPRESS}" ]]; then
  COMPRESS="f"
fi

if [[ -z "${CLEAN}" ]]; then
  CLEAN="f"
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  OUTPUT_DIR="mg-clust_output"
fi

if [[ -z "${NSLOTS}" ]]; then
  NSLOTS="12"
fi

if [[ -z "${NJOBS}" ]]; then
  NJOBS="3"
fi

if [[ -z "${MIN_CONTIG_LENGTH}" ]]; then
  MIN_CONTIG_LENGTH="250"
fi

if [[ -z "${MIN_OPU_OCCUP}" ]]; then
  MIN_OPU_OCCUP="2"
fi

if [[ -z "${MIN_ORF_LENGTH}" ]]; then
  MIN_ORF_LENGTH="60"
fi

if [[ -z "${OVERWRITE}" ]]; then
  OVERWRITE="f"
fi

if [[ -z "${RUN_MODULE_1}" ]]; then
  RUN_MODULE_1="t"
fi

if [[ -z "${RUN_MODULE_2}" ]]; then
  RUN_MODULE_2="t"
fi

if [[ -z "${RUN_MODULE_3}" ]]; then
  RUN_MODULE_3="t"
fi

if [[ -z "${RUN_MODULE_4}" ]]; then
  RUN_MODULE_4="t"
fi

if [[ -z "${SERVERS}" ]]; then
  SERVERS="lacalavera,oceania"
fi

if [[ -z "${TRAIN_FILE_NAME}" ]]; then
  TRAIN_FILE_NAME="illumina_1"
fi

if [[ -z "${THRES_RANGE}" ]]; then
  THRES_RANGE="0.7,0.9"
fi  

if [[ -z "${THRES_STEP}" ]]; then
  THRES_STEP="0.1"
fi  

###############################################################################
### 5. Check mandatory parameters
###############################################################################

NSLOTS_TOTAL=$(echo "${NJOBS}*${NSLOTS}" | bc)
N_SERVERS=$(echo "${SERVERS}" | awk -v FS="," '{print NF}')
MAX_JOBS=$(echo "${NJOBS}/${N_SERVERS}" | bc)
  
if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "${INPUT_DIR} is not a directory"
  exit 1
fi  

if [[ -z "${R1_SUFFIX}" ]]; then
  echo "missing suffix of R1 reads; Use --reads1_suffix"
  exit 1
fi  

if [[ -z "${R2_SUFFIX}" ]]; then
  echo "missing suffix of R2 reads; Use --reads2_suffix"
  exit 1
fi  

###############################################################################
### 6. Create output folder
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then

  if [[ "${OVERWRITE}" == "t" ]]; then
    
    if [[ "${RUN_MODULE_1}" == "t" && "${RUN_MODULE_2}" == "t" \
          && "${RUN_MODULE_3}" == "t" && "${RUN_MODULE_4}" == "t" ]]; then
    
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
    
    else 
    
      echo "Not all modules will be executed; it is not safe to overwrite"
      echo "To overwrite remove folder by hand"
      exit 0
    
    fi
    
  elif [[ "${OVERWRITE}" == "f" ]]; then
  
    if [[ "${RUN_MODULE_1}" == "t" && "${RUN_MODULE_2}" == "t" \
          && "${RUN_MODULE_3}" == "t" && "${RUN_MODULE_4}" == "t" ]]; then
    
      echo "${OUTPUT_DIR} already exists and all modules will be executed; use --overwrite t to overwrite"
      exit 0
      
    else
    
      echo "${OUTPUT_DIR} already exists; the data in ${OUTPUT_DIR} will be used to continue the analysis"
    
    fi
    
  fi
  
elif [[ ! -d "${OUTPUT_DIR}" ]]; then
  
   mkdir -p "${OUTPUT_DIR}"

  if [[ "$?" -ne "0" ]]; then
    echo "mkdir output directory ${OUTPUT_DIR} failed"
    exit 1
  fi  
fi

###############################################################################
### 7. Run module 1
###############################################################################

R1_LIST=$(ls "${INPUT_DIR}"/*"${R1_SUFFIX}")
R2_LIST=$(ls "${INPUT_DIR}"/*"${R2_SUFFIX}")
SAMPLE_NAME_LIST=$(echo "${R1_LIST}" | xargs -n 1 basename | sed "s/${R1_SUFFIX}//")

OUTPUT_DIR_1="${OUTPUT_DIR}/output-1"

if [[ "${RUN_MODULE_1}" == "t" ]]; then

  mkdir -p "${OUTPUT_DIR_1}"
  if [[ $? -ne "0" ]]; then
    echo "mkdir -p ${OUTPUT_DIR_1} failed"
    exit 1
  fi  

  if [[ -z "${ASSEM_DIR}" ]]; then
  
    env_parallel \
    -S "${SERVERS}" \
    --jobs "${NJOBS}" \
    --max-procs "${MAX_JOBS}" \
    --joblog "${LOGS_FILE}" \
    --link \
    "${MODULES}/mg-clust_module-1.bash \
    --reads1 {1} \
    --reads2 {2} \
    --assem_preset ${ASSEM_PRESET} \
    --nslots ${NSLOTS} \
    --output_dir ${OUTPUT_DIR_1}/{3} \
    --overwrite t \
    --sample_name {3} \
    --min_contig_length ${MIN_CONTIG_LENGTH} > ${OUTPUT_DIR_1}/module-1-{3}.log" \
    ::: "${R1_LIST}" ::: "${R2_LIST}" ::: "${SAMPLE_NAME_LIST}"
    
    if [[ "$?" -ne "0" ]]; then
      echo "mg-clust_module-1.bash failed"
      exit 1
    fi

  else 
  
  if [[ ! -d "${ASSEM_DIR}" ]]; then
    echo "${ASSEM_DIR} is not a real directory; please provide directory assembly as output-1/SAMPLE_NAME/SAMPLE_NAME.contigs.fa"
    exit 1
  fi
  
    env_parallel \
    -S "${SERVERS}" \
    --jobs "${NJOBS}" \
    --max-procs "${MAX_JOBS}" \
    --joblog "${LOGS_FILE}" \
    --link \
    "${MODULES}/mg-clust_module-1.bash \
    --reads1 {1} \
    --reads2 {2} \
    --assem_dir ${ASSEM_DIR} \
    --assem_preset ${ASSEM_PRESET} \
    --nslots ${NSLOTS} \
    --output_dir ${OUTPUT_DIR_1}/{3} \
    --overwrite t \
    --sample_name {3} \
    --min_contig_length ${MIN_CONTIG_LENGTH} > ${OUTPUT_DIR_1}/module-1-{3}.log" \
    ::: "${R1_LIST}" ::: "${R2_LIST}" ::: "${SAMPLE_NAME_LIST}"
  
    if [[ "$?" -ne "0" ]]; then
      echo "mg-clust_module-1.bash failed"
      exit 1
    fi
    
  fi

fi  

###############################################################################
### 8. Run module 2
###############################################################################

SAMPLE_NAME_LIST=$(ls -d "${OUTPUT_DIR_1}"/*/ | xargs -n1 basename | sed "s/\/$//")
OUTPUT_DIR_2="${OUTPUT_DIR}/output-2"

if [[ "${RUN_MODULE_2}" == "t" ]]; then

  mkdir -p "${OUTPUT_DIR_2}"
  if [[ $? -ne "0" ]]; then
    echo "mkdir -p ${OUTPUT_DIR_2} failed"
    exit 1
  fi  

  # check if there is an assembly file and reassign if not
  if [[ -z "${ASSEM_DIR}" ]]; then
    ASSEM_DIR_REASSIGNED="${OUTPUT_DIR_1}"
  else 
    ASSEM_DIR_REASSIGNED="${ASSEM_DIR}"
  fi 
    
  env_parallel \
  -S "${SERVERS}" \
  --jobs "${NJOBS}" \
  --max-procs "${MAX_JOBS}" \
  --joblog "${LOGS_FILE}" \
  "${MODULES}/mg-clust_module-2.bash \
  --assembly_file ${ASSEM_DIR_REASSIGNED}/{}/{}.contigs.fa \
  --bam_file ${OUTPUT_DIR_1}/{}/{}_sorted_markdup.bam \
  --nslots ${NSLOTS} \
  --output_dir ${OUTPUT_DIR_2}/{} \
  --overwrite t \
  --sample_name {} \
  --train_file_name ${TRAIN_FILE_NAME} > ${OUTPUT_DIR_2}/module-2-{}.log" \
  ::: "${SAMPLE_NAME_LIST}"

  if [[ "$?" -ne "0" ]]; then
    echo "mg-clust_module-2.bash failed"
    exit 1
  fi
  
fi  

###############################################################################
### 9. Run module 3
###############################################################################

OUTPUT_DIR_3="${OUTPUT_DIR}/output-3"

if [[ "${RUN_MODULE_3}" == "t" ]]; then

  mkdir "${OUTPUT_DIR_3}"
  if [[ $? -ne "0" ]]; then
    echo "mkdir ${OUTPUT_DIR_3} failed"
    exit 1
  fi 

  "${MODULES}/mg-clust_module-3.bash" \
  --input_dir "${OUTPUT_DIR_2}" \
  --min_orf_length "${MIN_ORF_LENGTH}" \
  --nslots "${NSLOTS}" \
  --output_dir "${OUTPUT_DIR_3}" > "${OUTPUT_DIR_3}/module-3.log"

  if [[ "$?" -ne "0" ]]; then
    echo "mg-clust_module-3.bash failed"
    exit 1
  fi  

fi

###############################################################################
### 10. Run module 4
###############################################################################

OUTPUT_DIR_4="${OUTPUT_DIR}/output-4"

if [[ "${RUN_MODULE_4}" == "t" ]]; then

  if [[ ! -d "${OUTPUT_DIR_4}" ]]; then
  
    mkdir "${OUTPUT_DIR_4}"
    
    if [[ "$?" -ne "0" ]]; then
      echo "mkdir  ${OUTPUT_DIR_4} failed"
      exit 1
    fi
  
  else 
    
    echo "${OUTPUT_DIR_4} already exists; this folder will be used to continue the analysis"
  
  fi
    
  THRES_MIN=$(echo ${THRES_RANGE} | cut -f1 -d",")
  THRES_MAX=$(echo ${THRES_RANGE} | cut -f2 -d",")
  THRES_SEQ=$(LC_ALL=C seq "${THRES_MIN}" "${THRES_STEP}" "${THRES_MAX}")
               
  ORFS_DB="${OUTPUT_DIR_3}/orfs_filtered_db"
  COV_TABLE="${OUTPUT_DIR_3}/orfs_cov.tsv"
  
  env_parallel \
  -S "${SERVERS}" \
  --jobs "${NJOBS}" \
  --max-procs "${MAX_JOBS}" \
  --joblog "${LOGS_FILE}" \
  "${MODULES}/mg-clust_module-4.bash \
  --orfs_db ${ORFS_DB} \
  --output_dir ${OUTPUT_DIR_4} \
  --cov_table "${COV_TABLE}" \
  --nslots "${NSLOTS}" \
  --min_opu_occup ${MIN_OPU_OCCUP} \
  --clust_thres {} > ${OUTPUT_DIR_4}/module-4-{}.log" \
  ::: $(echo "${THRES_SEQ}")

  if [[ "$?" -ne "0" ]]; then
    echo "mg-clust_module-4.bash failed"
    exit 1
  fi  

fi

###############################################################################
### 11. Clean
###############################################################################

if [[ "${CLEAN}" == "t" && "${RUN_MODULE_1}" == "t" ]]; then
  rm -r "${OUTPUT_DIR_1}" 
  if [[ $? -ne "0" ]]; then
    echo "rm -r ${OUTPUT_DIR_1} (intermediate data) failed"
    exit 1
  fi
fi

if [[ "${CLEAN}" == "t" && "${RUN_MODULE_2}" == "t" ]]; then
  rm -r "${OUTPUT_DIR_2}" 
  if [[ $? -ne "0" ]]; then
    echo "rm -r ${OUTPUT_DIR_2} (intermediate data) failed"
    exit 1
  fi
fi

if [[ "${CLEAN}" == "t" && "${RUN_MODULE_3}" == "t" ]]; then
  rm -r "${OUTPUT_DIR_3}" 
  if [[ $? -ne "0" ]]; then
    echo "rm -r ${OUTPUT_DIR_3} (intermediate data) failed"
    exit 1
  fi
fi


if [[ "${CLEAN}" == "t" && "${RUN_MODULE_4}" == "t" ]]; then
  rm -r "${OUTPUT_DIR_4}" 
  if [[ $? -ne "0" ]]; then
    echo "rm -r ${OUTPUT_DIR_4} (intermediate data) failed"
    exit 1
  fi
fi

###############################################################################
### 12. Compress
###############################################################################

if [[ "${COMPRESS}" == "t" ]]; then

  if [[ "${CLEAN}" == "f" && "${RUN_MODULE_1}" == "t" ]]; then
  
    tar czfv "${OUTPUT_DIR}/output-1.tar.gz" \
    --directory $(dirname "${OUTPUT_DIR_1}") ./output-1 \
    --remove-files 
    
    if [[ $? -ne "0" ]]; then
      echo "tar czfv output-1 failed"
      exit 1
    fi
    
  fi
  
  if [[ "${CLEAN}" == "f" && "${RUN_MODULE_2}" == "t" ]]; then
  
    tar czfv "${OUTPUT_DIR}/output-2.tar.gz" \
    --directory $(dirname "${OUTPUT_DIR_2}") ./output-2 \
    --remove-files 
  
    if [[ $? -ne "0" ]]; then
      echo "tar czfv output-2 failed"
      exit 1
    fi
  
  fi
  
  if [[ "${CLEAN}" == "f" && "${RUN_MODULE_3}" == "t" ]]; then
  
    tar czfv "${OUTPUT_DIR}/output-3.tar.gz" \
    --directory $(dirname "${OUTPUT_DIR_3}") ./output-3 \
    --remove-files 
  
    if [[ $? -ne "0" ]]; then
      echo "tar czfv output-3 failed"
      exit 1
    fi
  
  fi
  
  if [[ "${CLEAN}" == "f" && "${RUN_MODULE_4}" == "t" ]]; then
  
    tar czfv "${OUTPUT_DIR}/output-4.tar.gz" \
    --directory $(dirname "${OUTPUT_DIR_4}") ./output-4 \
    --remove-files 
  
    if [[ $? -ne "0" ]]; then
      echo "tar czfv output-4 failed"
      exit 1
    fi
  
  fi
  
  if [[ $(ls "${OUTPUT_DIR}/orfs_clust2abund_id"*".tsv" 2> /dev/null | wc -l ) -gt 0 ]]; then
  
    env_parallel \
    -S "${SERVERS}" \
    -j "${NJOBS}" \
    --joblog "${LOGS_FILE}" \
    "pigz --processes ${NSLOTS} {}" \
    ::: $(ls "${OUTPUT_DIR}/orfs_clust2abund_id"*".tsv")

    if [[ $? -ne "0" ]]; then
      echo "pigz --processes ${NSLOTS} ${OUTPUT_DIR}/orfs_clust2abund_id*.tsv failed"
      exit 1
    fi
    
  fi
  
fi

###############################################################################
### 13. Exit
###############################################################################

echo "mg-clust exited successfully"
exit 0


