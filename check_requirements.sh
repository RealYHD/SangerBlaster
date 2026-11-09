#!/bin/bash
# check_requirements.sh - Deterministic pre-flight requirement check

set -e  # Exit on any error

echo "=== GenoSanger Pipeline Requirement Check ==="
echo ""

# Check container runtime
check_container() {
    if command -v docker &> /dev/null; then
        echo "✓ Docker is available"
        docker --version | head -1
        return 0
    elif command -v singularity &> /dev/null; then
        echo "✓ Singularity is available"
        singularity version 2>&1 | head -1
        return 0
    else
        echo "✗ ERROR: Neither Docker nor Singularity is installed!"
        echo ""
        echo "To install:"
        echo "  Docker:   https://docs.docker.com/engine/install/"
        echo "  Singularity: apt-get install singularity-container"
        exit 1
    fi
}

# Check download tools
check_download_tools() {
    if command -v wget &> /dev/null || command -v curl &> /dev/null; then
        echo "✓ wget or curl is available"
        return 0
    else
        echo "✗ ERROR: Neither wget nor curl is available!"
        echo ""
        echo "To install:"
        echo "  apt-get install wget curl"
        exit 1
    fi
}

# Check Nextflow version (optional but recommended)
check_nextflow() {
    if command -v nextflow &> /dev/null; then
        version=$(nextflow -version | head -1)
        echo "✓ Nextflow is available: $version"
        
        # Check if version is too old
        current_version=$(nextflow -version | grep -oP '(\d+\.\d+\.\d+)')
        min_version="24.04.6"  # Recommended minimum
        
        if [ "$(printf '%s\n' "$min_version" "$current_version" | sort -V | head -n1)" = "$min_version" ]; then
            echo "✓ Nextflow version is acceptable ($current_version ≥ $min_version)"
        else
            echo "⚠ Warning: Nextflow $current_version may be too old."
            echo "  Recommended: $min_version or newer"
            echo "  Download: https://get.nextflow.io/nf/"
        fi
        return 0
    else
        echo "✗ ERROR: Nextflow is not installed!"
        echo ""
        echo "To install:"
        echo "  curl -s https://get.nextflow.io/nf/ | bash"
        echo "  or: wget https://github.com/nextflow-io/nextflow/releases/download/v26.10.5/nextflow -O nextflow && chmod +x nextflow"
        exit 1
    fi
}

# Main check
echo "Checking requirements..."
check_container
check_download_tools
check_nextflow

echo ""
echo "=== All requirements met! Pipeline can proceed ==="
exit 0
