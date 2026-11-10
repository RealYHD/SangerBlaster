#!/bin/bash -e
# Build GenoSanger container with Tracy pre-installed
# Run: ./build_genosanger_container.sh

BASE_IMAGE="quay.io/nf-core/ubuntu:20.04"
CONTAINER_NAME="${1:-genosanger_ubuntu_20.04.img}"
TRACY_VERSION="0.9.1"
TRACY_URL="https://github.com/gear-genomics/tracy/releases/download/v${TRACY_VERSION}/tracy-v${TRACY_VERSION}-linux-amd64"

echo "=== Building GenoSanger container with Tracy ==="
echo "Base image: ${BASE_IMAGE}"
echo "Output: ${CONTAINER_NAME}"

# Pull base Singularity image
echo "Pulling base image..."
singularity pull --name "${CONTAINER_NAME}" "${BASE_IMAGE}"

# Check if we already have a build directory
BUILD_DIR=$(mktemp -d -t genosanger-build-XXXXXXXXXX)
trap "rm -rf ${BUILD_DIR}" EXIT

echo "Building container with Tracy pre-installed..."

# Create definition file
cat > "${BUILD_DIR}/genosanger.def" <<'EOF'
Bootstrap: docker
From: quay.io/nf-core/ubuntu:20.04

%post
    # Update apt and install curl (pre-installed in Ubuntu, but ensure it)
    apt-get update -qq
    apt-get install -y -qq curl
    
    # Download Tracy binary
    curl -L -o /usr/local/bin/tracy "${TRACY_URL}"
    chmod +x /usr/local/bin/tracy
    
    # Verify Tracy installation
    tracy --version

%environment
    export TRACY_PATH=/usr/local/bin/tracy

%runscript
    exec /bin/bash
EOF

# Build Singularity image from definition
singularity build "${CONTAINER_NAME}" "${BUILD_DIR}/genosanger.def"

echo "=== Container built successfully ==="
echo "Location: $(pwd)/${CONTAINER_NAME}"
echo "Test it:"
echo "  singularity run ${CONTAINER_NAME} tracy --version"
