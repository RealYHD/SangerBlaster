#!/bin/bash
# check_requirements.sh - Deterministic pre-flight requirement check

set -e  # Exit on any error

echo "=== GenoSanger Pipeline Requirement Check ==="
echo ""

# Check Conda/Mamba environment (recommended)
check_conda() {
    if command -v conda &> /dev/null || command -v mamba &> /dev/null; then
        echo "✓ Conda/Mamba is available"
        
        # Check for tracy package
        if mamba list 2>/dev/null | grep -q tracy || conda list 2>/dev/null | grep -q tracy; then
            echo "✓ Tracy is installed in environment"
            
            # Find which environment has tracy
            if command -v mamba &> /dev/null; then
                env=$(mamba list -f regex 'tracy' 2>/dev/null | head -1 | awk '{print $2}')
                echo "  Environment: $env"
            else
                env=$(conda list | grep tracy)
                echo "  Package: $env"
            fi
            
            # Verify tracy command works
            if command -v tracy &> /dev/null; then
                tracy --version | head -1
                echo "✓ Tracy version OK"
            else
                echo "⚠ Warning: tracy package found but command not available"
                echo "  Try: source activate genosanger"
            fi
            
            return 0
        else
            echo "✗ ERROR: Tracy package not installed in Conda environment!"
            echo ""
            echo "Install with:"
            echo "  mamba create -y -n genosanger -c conda-forge -c bioconda tracy blast"
            echo "  source activate genosanger"
            exit 1
        fi
    else
        echo "✗ ERROR: Conda/Mamba is not installed!"
        echo ""
        echo "To install:"
        echo "  Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        echo "  Mamba: conda install -y mamba -c conda-forge"
        exit 1
    fi
}

# Fallback: Check container runtime (if Conda not available)
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
        return 1
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

if check_conda; then
    # Conda/tracy is available - this is the recommended path
    check_nextflow
elif check_container; then
    # Container fallback - but warn that this is less preferred
    echo "⚠ Warning: Using container runtime instead of Conda"
    check_nextflow
else
    echo "✗ ERROR: Neither Conda/tracy nor container runtime available!"
    exit 1
fi

echo ""
echo "=== All requirements met! Pipeline can proceed ==="
exit 0
