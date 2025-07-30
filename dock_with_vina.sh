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
PDBID="$4"  # PDB ID
AUR="$5"  # Aurora version
EXP="$6"  # Experiment name


if [[ "$PDBID" == "4af3" ]] || [[ "$AUR" == "B" ]]; then
    RECEPTOR="$SCRIPT_DIR/4af3/4af3_A_rec_reduce_noflip.pdb"
    RECEPTOR_PREFIX="4af3"
    CENTER="21 -21 12"
elif [[ "$PDBID" == "4ceg" ]] || [[ "$AUR" == "A" ]]; then
    RECEPTOR="$SCRIPT_DIR/4ceg/4ceg_1_protein.pdb"
    RECEPTOR_PREFIX="4ceg_receptor"
    CENTER="10 20 5"
# elif [[ "$AURORA" == "C" ]]; then
#     RECEPTOR="$SCRIPT_DIR/1iep/1iep.pdb"
#     RECEPTOR_PREFIX="1iep"
#     CENTER="15.190 53.903 16.917"
else
    echo "Error: --pdbid argument must be '4af3' or '4ceg' (or aurora A/B)"
    echo "Received PDBID='$PDBID', AURORA='$AUR'"
    exit 1
fi

