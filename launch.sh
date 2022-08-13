#!/bin/bash
#
# Copyright (c) 2022, NVIDIA CORPORATION.
# SPDX-License-Identifier: Apache-2.0

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

LOCAL_ENV=.env


function config() {
    local write_env=$1

    DOCKER_IMAGE="flow:latest"
    PROJECT_PATH=${PROJECT_PATH:=$(pwd)}
    DATA_PATH=${DATA_PATH:=$(pwd)/data}

    if [ $write_env -eq 1 ]; then
        echo DOCKER_IMAGE=${DOCKER_IMAGE} > $LOCAL_ENV
        echo PROJECT_PATH=${PROJECT_PATH} >> $LOCAL_ENV
        echo DATA_PATH=${DATA_PATH} >> $LOCAL_ENV
    fi
}


build() {
    echo -e "Building ${DOCKER_IMAGE}..."
    DOCKER_BUILD_CMD="docker build \
        -t ${DOCKER_IMAGE} \
        -f docker/Dockerfile.flow"

    while [[ $# -gt 0 ]]; do
        case $1 in
            -b|--base-image)
                BASE_IMAGE=$2
                shift
                shift
                ;;
            -c|--clean)
                DOCKER_BUILD_CMD="${DOCKER_BUILD_CMD} --no-cache"
                shift
                ;;
            *)
                echo "Unknown option $1."
                exit 1
                ;;
        esac
    done

    if [ ! -z "${BASE_IMAGE}" ];
    then
        DOCKER_BUILD_CMD="${DOCKER_BUILD_CMD} --build-arg BASE_IMAGE=${BASE_IMAGE}"
    fi

    $DOCKER_BUILD_CMD .
}


dev() {
    local CMD='bash'

    while [[ $# -gt 0 ]]; do
        case $1 in
            -a|--additional-args)
                DOCKER_CMD="${DOCKER_CMD} $2"
                shift
                shift
                ;;
            -i|--image)
                DEV_IMG="$2"
                shift
                shift
                ;;
            -d|--demon)
                DOCKER_CMD="${DOCKER_CMD} -d"
                shift
                ;;
            -c|--cmd)
                shift
                CMD="$@"
                break
                ;;
            *)
                echo "Unknown option '$1'.
Available options are -a(--additional-args), -i(--image), -d(--demon) and -c(--cmd).
Please always ensure -c is the last option."
                exit 1
                ;;
        esac
    done

    set -x
    $DOCKER_CMD --rm ${DEV_IMG} ${CMD}
}


if [ -e ./${LOCAL_ENV} ]
then
    echo -e "sourcing environment from ./${LOCAL_ENV}"
    . ./${LOCAL_ENV}
    config 0
else
    echo -e "${YELLOW}Writing deafults to ${LOCAL_ENV}${RESET}"
    config 1
fi

DOCKER_CMD="docker run --gpus all --network=host\
    -v ${DATA_PATH}:/data\
    -v ${PROJECT_PATH}:/workspace\
    -v ${CHEMBLE_DB}:/data/chembl.db\
    -e PYTHONPATH=/workspace\
    -it"


case $1 in
    config)
        config 1
        ;;
    build)
        "$@"
        ;;
    dev)
        "$@"
        ;;
    *)
        usage
        ;;
esac
