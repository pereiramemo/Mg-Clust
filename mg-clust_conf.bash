#dirs
#BIOINF_BIN="/home/bioinf/bin"
BIOINF_BIN="/nfs/bin"
HOME="/home/epereira"
REPOS="${HOME}/workspace/repositories"
#MODULES="${REPOS}/tools/metagenomic_pipelines/metaclust/modules"
MODULES="${BIOINF_BIN}/metaclust/modules"

# tools
megahit_version="1.2.9-Linux-x86_64-static"
megahit="${BIOINF_BIN}/megahit/MEGAHIT-${megahit_version}/bin/megahit"

bwa_version="0.7.17-r1198-dirty"
bwa="${BIOINF_BIN}/bwa/bwa-${bwa_version}/bwa"

samtools_version="1.9"
samtools="${BIOINF_BIN}/samtools/samtools-${samtools_version}/samtools"

picard="${BIOINF_BIN}/picard/picard.jar"

# fraggenescan_version="-1.31"
# fraggenescan="${BIOINF_BIN}/fraggenescan/FragGeneScan${fraggenescan_version}/FragGeneScan"

fraggenescan_version="Rs"
fraggenescan="${BIOINF_BIN}/fraggenescan/FragGeneScan${fraggenescan_version}/bin/FragGeneScan${fraggenescan_version}"
FRAGGENESCAN_TRAIN_FILE_DIR="${BIOINF_BIN}/fraggenescan/FragGeneScan${fraggenescan_version}/FragGeneScan${fraggenescan_version}/train"


bedtools_version="v2.30.0"          
bedtools="${BIOINF_BIN}/bedtools/bedtools_${bedtools_version}/bin/bedtools"

bbduk_version="38.79"
bbduk="${BIOINF_BIN}/bbmap/bbmap-${bbduk_version}/bbduk.sh"

mmseqs="${BIOINF_BIN}/mmseqs2/mmseqs/bin/mmseqs"

source "/usr/bin/env_parallel.bash"
