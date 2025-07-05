# config.sh - Configuration file for run_pipeline.sh

# === Set script directory ===
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script directory: $SCRIPT_DIR"

# === Configuration ===
# Default values
RECEPTOR=""
RECEPTOR_PREFIX=""
CENTER=""
SIZE="30 30 30"
EXHAUST="32"

EPOCH="$1"  # Epoch number
MOLS="$2"  # Number of molecules
BS="$3"    # Batch size
PDBID="$4"  # Aurora version


if [[ "$PDBID" == "4af3" ]]; then
    RECEPTOR="$SCRIPT_DIR/4af3/4af3.pdb"
    RECEPTOR_PREFIX="4af3"
    CENTER="21 -21 12"
elif [[ "$PDBID" == "A" ]]; then
    RECEPTOR="$SCRIPT_DIR/4ceg/4ceg_1_protein.pdb"
    RECEPTOR_PREFIX="4ceg_receptor"
    CENTER="10 20 5"
# elif [[ "$AURORA" == "C" ]]; then
#     RECEPTOR="$SCRIPT_DIR/1iep/1iep.pdb"
#     RECEPTOR_PREFIX="1iep"
#     CENTER="15.190 53.903 16.917"
else
    echo "Error: --pdbid argument must be A or B"
    exit 1
fi

