#!/usr/bin/env bash
#
# clean_nextflow_cache.sh
#
# Cleanup script to remove Nextflow generated work directories and cache
# Run this before starting a fresh pipeline run or to free up disk space.
#

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

# Default directories to clean
NEXTFLOW_CACHE="${HOME}/.nextflow"
CONDA_CACHE="${NEXTFLOW_CACHE}/conda/cache"
WORK_DIR="${PROJECT_ROOT}/work"

echo -e "${YELLOW}Nextflow Cache Cleanup Script${NC}"
echo "=================================="
echo ""

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --dry-run    Show what would be deleted without actually deleting"
    echo "  --force      Skip confirmation prompts"
    echo "  --help       Show this help message"
    echo ""
    echo "Default directories to clean:"
    echo "  - ~/.nextflow (Nextflow system cache)"
    echo "  - ~/.nextflow/conda/cache (Conda environment cache)"
    echo "  - ./work (Nextflow work directory in project root)"
}

# Parse arguments
DRY_RUN=false
FORCE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --help)
            show_usage
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            show_usage
            exit 1
            ;;
    esac
done

# Function to list and remove directories
cleanup_dir() {
    local dir="$1"
    local name="$2"
    
    if [[ -d "$dir" ]]; then
        if [[ "$DRY_RUN" == true ]]; then
            echo -e "${YELLOW}[DRY RUN]${NC} Would remove: $name ($dir)"
            # List contents for dry run
            if [[ "$name" == "Nextflow system cache" ]]; then
                echo "  Contents:"
                find "$dir" -maxdepth 2 -type f | head -20 | sed 's/^/    /'
            fi
        else
            echo "Removing: $name ($dir)"
            rm -rf "$dir"
            echo -e "${GREEN}✓ Removed${NC}"
        fi
    else
        echo "Skipped: $name ($dir) - does not exist"
    fi
}

echo ""
echo "Cleaning up Nextflow directories..."
echo ""

# Clean each directory
cleanup_dir "$NEXTFLOW_CACHE" "Nextflow system cache"
cleanup_dir "$CONDA_CACHE" "Conda environment cache"
cleanup_dir "$WORK_DIR" "Project work directory"

echo ""
echo -e "${YELLOW}Note: Conda cache cleanup may take a while...${NC}"
echo ""

if [[ "$DRY_RUN" == true ]]; then
    echo -e "${YELLOW}[DRY RUN] Complete. No files were actually deleted.${NC}"
else
    if [[ "$FORCE" == false ]]; then
        read -p "Press Enter to continue..." || true
    fi
    echo -e "${GREEN}Cleanup completed successfully!${NC}"
fi
