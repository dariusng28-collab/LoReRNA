#!/usr/bin/env bash
# =============================================================================
# build_containers.sh — build and push LoReRNA Docker images
#
# Run this on a LOCAL MACHINE that has Docker and is logged in to Docker Hub.
# Do NOT run on HPC nodes — Docker is typically unavailable there.
# On HPC, Nextflow pulls the images from Docker Hub and converts to SIF
# automatically using Singularity.
#
# Usage:
#   cd /path/to/lorerna           # must run from repo root
#   bash containers/build_containers.sh [--push] [--version 1.1.0]
#
# Flags:
#   --push           push to Docker Hub after building (default: build only)
#   --version TAG    image tag to build (default: reads from nextflow.config)
#   --miser-only     rebuild only the MisER image
#   --lorerna-only   rebuild only the lorerna image
#   --platform       multi-arch build (default: linux/amd64,linux/arm64)
#                    requires: docker buildx create --use
#
# Singularity note:
#   After pushing, users on HPC do NOT need to manually pull.
#   Nextflow's singularity profile pulls automatically into singularity.cacheDir
#   the first time each image is used. To force a manual pull:
#
#   singularity pull \
#     /path/to/singularity/cache/dariusng28-lorerna-1.1.0.sif \
#     docker://dariusng28/lorerna:1.1.0
#
#   singularity pull \
#     /path/to/singularity/cache/dariusng28-miser-1.1.0.sif \
#     docker://dariusng28/miser:1.1.0
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ── Parse flags ───────────────────────────────────────────────────────────────
PUSH=false
VERSION=""
BUILD_LORERNA=true
BUILD_MISER=true
PLATFORM="linux/amd64"   # change to "linux/amd64,linux/arm64" for multi-arch

while [[ $# -gt 0 ]]; do
    case "$1" in
        --push)         PUSH=true;          shift ;;
        --version)      VERSION="$2";       shift 2 ;;
        --miser-only)   BUILD_LORERNA=false; shift ;;
        --lorerna-only) BUILD_MISER=false;   shift ;;
        --platform)     PLATFORM="$2";      shift 2 ;;
        *) echo "Unknown flag: $1"; exit 1 ;;
    esac
done

# ── Resolve version from container_lorerna tag in nextflow.config ────────────
if [ -z "${VERSION}" ]; then
    VERSION=$(grep "container_lorerna" "${REPO_ROOT}/nextflow.config" \
              | grep -v '^\s*//' | grep -oP "lorerna:\K[\d]+\.[\d]+\.[\d]+")
    if [ -z "${VERSION}" ]; then
        echo "ERROR: Could not parse container version from nextflow.config."
        echo "       Expected a line like: container_lorerna = 'docker://user/lorerna:1.1.0'"
        echo "       Set it explicitly: --version 1.1.0"
        exit 1
    fi
fi

DOCKER_USER="dariusng28"
LORERNA_TAG="${DOCKER_USER}/lorerna:${VERSION}"
MISER_TAG="${DOCKER_USER}/miser:${VERSION}"

echo "============================================"
echo " LoReRNA container build"
echo " Version     : ${VERSION}"
echo " lorerna tag : ${LORERNA_TAG}"
echo " miser tag   : ${MISER_TAG}"
echo " Platform    : ${PLATFORM}"
echo " Push        : ${PUSH}"
echo " Repo root   : ${REPO_ROOT}"
echo "============================================"

# ── Build lorerna (samtools, isoquant, oarfish, R, MultiQC) ──────────────────
if ${BUILD_LORERNA}; then
    echo
    echo "[1/2] Building lorerna image: ${LORERNA_TAG}"
    echo "      context : ${REPO_ROOT}"
    echo "      file    : ${REPO_ROOT}/containers/Dockerfile"

    if ${PUSH}; then
        docker buildx build \
            --platform "${PLATFORM}" \
            --file     "${REPO_ROOT}/containers/Dockerfile" \
            --tag      "${LORERNA_TAG}" \
            --push \
            "${REPO_ROOT}"
    else
        docker build \
            --platform "${PLATFORM}" \
            --file     "${REPO_ROOT}/containers/Dockerfile" \
            --tag      "${LORERNA_TAG}" \
            "${REPO_ROOT}"
    fi

    echo "lorerna build complete: ${LORERNA_TAG}"
fi

# ── Build miser (MisER, pysam, samtools=1.23.1) ───────────────────────────────
if ${BUILD_MISER}; then
    echo
    echo "[2/2] Building miser image: ${MISER_TAG}"
    echo "      context : ${REPO_ROOT}"
    echo "      file    : ${REPO_ROOT}/containers/Dockerfile.miser"

    if ${PUSH}; then
        docker buildx build \
            --platform "${PLATFORM}" \
            --file     "${REPO_ROOT}/containers/Dockerfile.miser" \
            --tag      "${MISER_TAG}" \
            --push \
            "${REPO_ROOT}"
    else
        docker build \
            --platform "${PLATFORM}" \
            --file     "${REPO_ROOT}/containers/Dockerfile.miser" \
            --tag      "${MISER_TAG}" \
            "${REPO_ROOT}"
    fi

    echo "miser build complete: ${MISER_TAG}"
fi

# ── Post-build validation ─────────────────────────────────────────────────────
if ! ${PUSH}; then
    echo
    echo "============================================"
    echo " Build complete (not pushed — run with --push to publish)"
    echo " Quick sanity check:"
    echo "============================================"

    if ${BUILD_LORERNA}; then
        echo
        echo " lorerna:"
        docker run --rm "${LORERNA_TAG}" samtools  --version | head -1
        docker run --rm "${LORERNA_TAG}" isoquant --version 2>&1 | head -1 || true
        docker run --rm "${LORERNA_TAG}" oarfish --version 2>&1 | head -1 || true
        docker run --rm "${LORERNA_TAG}" Rscript --version 2>&1 | head -1
        docker run --rm "${LORERNA_TAG}" multiqc --version
    fi

    if ${BUILD_MISER}; then
        echo
        echo " miser:"
        docker run --rm "${MISER_TAG}" samtools --version | head -1
        docker run --rm "${MISER_TAG}" python3 -c "import MisER; print('MisER OK')" 2>/dev/null || \
            docker run --rm "${MISER_TAG}" MisER --version 2>&1 | head -1 || true
    fi

    echo
    echo " To push to Docker Hub:"
    echo "   bash containers/build_containers.sh --push --version ${VERSION}"
    echo
    echo " To manually pull on HPC after pushing:"
    if ${BUILD_LORERNA}; then
        echo "   singularity pull /path/to/singularity/cache/dariusng28-lorerna-${VERSION}.sif \\"
        echo "     docker://${LORERNA_TAG}"
    fi
    if ${BUILD_MISER}; then
        echo "   singularity pull /path/to/singularity/cache/dariusng28-miser-${VERSION}.sif \\"
        echo "     docker://${MISER_TAG}"
    fi
else
    echo
    echo "============================================"
    echo " Pushed to Docker Hub:"
    ${BUILD_LORERNA} && echo "   https://hub.docker.com/r/${DOCKER_USER}/lorerna/tags"
    ${BUILD_MISER}   && echo "   https://hub.docker.com/r/${DOCKER_USER}/miser/tags"
    echo
    echo " Update nextflow.config to reference the new tags:"
    echo "   container_lorerna = 'docker://${LORERNA_TAG}'"
    echo "   container_miser   = 'docker://${MISER_TAG}'"
    echo "============================================"
fi
