#!/bin/bash -l
#SBATCH --account ent_aiapps_omics
#SBATCH --partition batch_dgx1_m3
#SBATCH --nv-meta ml-model:virtual_screening,dcgm_opt_out:no
#SBATCH --nodes 1
#SBATCH --gpus-per-node 8
#SBATCH --ntasks-per-node 1
#SBATCH --time=8:00:00
#SBATCH --output run.log
#SBATCH --mem=0                 # all mem avail
#SBATCH --mail-type=ALL         # only send email on failure
#SBATCH --overcommit            # Needed for pytorch
#SBATCH --exclusive             # exclusive node access

### sbatch --partition <= removed. Do not know exact value
### sbatch --mail-user <userid> -J 1 --output train_1.log scripts/run_training_dp.sh
### sbatch --mail-user <userid> -J 1 --output train_1.log --nodes 4 --gpus-per-node 8 --ntasks 32 --ntasks-per-node 8  scripts/run_training_dp.sh


execute() {
    local RUN_CONTAINER=$1
    local CORRELATION_ID=$2
    local TASK_NAME=$3

    echo "Executing task '$TASK_NAME' belonging to job '$CORRELATION_ID' on container '${PAYLOAD_CONTAINER}'..."
    local JOB_DIR=/gpfs/fs1/projects/ent_aiapps/users/${USER}/data/${CORRELATION_ID}
    local TASK_DIR=${JOB_DIR}/${TASK_NAME}

    if [ ! -d ${JOB_DIR} ]; then
        echo "Job directory ${JOB_DIR} does not exist. Please create it before execute."
        return
    fi

    # Container Paths for mapping
    local DATA_PATH=/data
    local CONFIG_PATH=/config
    local OUTPUT_PATH=/output

    # Corresponding host paths. Mounting JOB_DIR instead of TASK_DIR to share
    # data between different tasks in a pipeline/job
    local HOST_DATA_PATH=${JOB_DIR}${DATA_PATH}
    local HOST_CONFIG_PATH=${JOB_DIR}${CONFIG_PATH}
    local HOST_OUTPUT_PATH=${JOB_DIR}${OUTPUT_PATH}

    local MOUNTS="${HOST_DATA_PATH}:${DATA_PATH},${HOST_CONFIG_PATH}:${CONFIG_PATH},${HOST_OUTPUT_PATH}:${OUTPUT_PATH}"

    # Create required directories and start the task
    echo "$TASK_NAME,SUBMITTING" > ${JOB_DIR}/status
    srun \
        --mpi=pmix \
        --nodes ${SLURM_JOB_NUM_NODES} \
        --ntasks ${SLURM_NTASKS} \
        --ntasks-per-node ${SLURM_NTASKS_PER_NODE} \
        --gpus-per-node ${SLURM_GPUS_PER_NODE} \
        --container-image ${RUN_CONTAINER} \
        --container-mounts ${MOUNTS} \
        --container-workdir ${WORKDIR} \
        python --version
}


usage() {
    cat <<EOF
USAGE: launch.sh
SLRUM launch job utility
----------------------------------------
launch.sh [options]

options:
    --id or -i
        Pipeline instance id. This will be a unique identifier for each job.
    --task or -t
        Task name in the pipeline
EOF
}


# PAYLOAD_CONTAINER="nvcr.io#nvidia/clara/megamolbart:0.1.2"
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --id | -i)
            CORRELATION_ID=$2
            shift
            shift
            ;;
        --task | -t)
            TASK_NAME=$2
            shift
            shift
            ;;
        --container | -c)
            PAYLOAD_CONTAINER=$2
            shift
            shift
            ;;
        *)
            usage
            exit
            ;;
    esac
done

if [[ ! -z ${PAYLOAD_CONTAINER} && ! -z ${CORRELATION_ID} && ! -z ${TASK_NAME} ]]; then
    execute "${PAYLOAD_CONTAINER}" "${CORRELATION_ID}" "${TASK_NAME}"
else
    echo "Please pass values for --id and --task options."
    exit 1
fi