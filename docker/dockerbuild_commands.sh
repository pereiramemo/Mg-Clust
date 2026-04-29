SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-1 -t ghcr.io/epereira/mg-clust/module-1:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-2 -t ghcr.io/epereira/mg-clust/module-2:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-3-4 -t ghcr.io/epereira/mg-clust/module-3-4:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-5 -t ghcr.io/epereira/mg-clust/module-5:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-6 -t ghcr.io/epereira/mg-clust/module-6:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-7 -t ghcr.io/epereira/mg-clust/module-7:latest .