#!/usr/bin/env bash
# Pre-run check for Conda/Tracy availability

echo "Checking Conda and Tracy installation..."

# Check for conda/mamba
if command -v conda &> /dev/null; then
    echo "✓ Conda found: $(conda --version)"
    CONDA_AVAILABLE=true
elif command -v mamba &> /dev/null; then
    echo "✓ Mamba found: $(mamba --version)"
    CONDA_AVAILABLE=true
else
    echo "✗ ERROR: Neither conda nor mamba found in PATH"
    echo "Please install Miniconda or Mamba:"
    echo "  wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3"
    echo "  export PATH=$HOME/miniconda3/bin:\$PATH"
    exit 1
fi

# Check Bioconda channel availability
echo "Checking Bioconda channel..."
conda config --add --force --channel bioconda 2>/dev/null || true
if conda search -c bioconda tracy &> /dev/null; then
    echo "✓ Tracy package available in Bioconda"
else
    echo "✗ ERROR: Could not find 'tracy' in Bioconda channel"
    echo "Try manually: conda search -c bioconda tracy"
    exit 1
fi

# Check Nextflow version
echo "Checking Nextflow version..."
if nextflow -version 2>&1 | grep -q "26."; then
    echo "✓ Nextflow 26.x found: $(nextflow -version | head -1)"
else
    echo "⚠ Nextflow version may be incompatible (expected 26.x)"
fi

echo ""
echo "All checks passed! Ready to run pipeline."
