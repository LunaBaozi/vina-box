#!/bin/bash
set -e  # Exit on error
set -u  # Error on undefined variables

# === Source configuration file (variables only, not environment) ===
# Don't source dock_with_vina.sh as it might interfere with conda environment
# source "$(dirname "${BASH_SOURCE[0]}")/dock_with_vina.sh"

# === Set variables manually (from dock_with_vina.sh) ===
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script directory: $SCRIPT_DIR"

# Set variables that would normally come from dock_with_vina.sh
EPOCH="$1"
MOLS="$2" 
BS="$3"
PDBID="$4"
AUR="$5"
EXP="$6"

# Set receptor configuration based on PDBID/Aurora
if [[ "$PDBID" == "4af3" ]] || [[ "$AUR" == "B" ]]; then
    RECEPTOR="$SCRIPT_DIR/4af3/4af3.pdb"   #4af3_A_rec_reduce_noflip.pdb"
    RECEPTOR_PREFIX="4af3"
    CENTER="21 -21 12"
elif [[ "$PDBID" == "4ceg" ]] || [[ "$AUR" == "A" ]]; then
    RECEPTOR="$SCRIPT_DIR/4ceg/4ceg_1_protein.pdb"
    RECEPTOR_PREFIX="4ceg_receptor"
    CENTER="10 20 5"
else
    echo "Error: --pdbid argument must be '4af3' or '4ceg' (or aurora A/B)"
    echo "Received PDBID='$PDBID', AURORA='$AUR'"
    exit 1
fi

SIZE="40 40 40"
EXHAUST="8"

BASE_DIR="$SCRIPT_DIR/$RECEPTOR_PREFIX/experiment_${EXP}_${EPOCH}_${MOLS}_${BS}_${PDBID}"
LIGAND_DIR="$BASE_DIR/ligands"
PREPARED_LIG_DIR="$BASE_DIR/prepared_ligands"
PREPARED_REC_DIR="$SCRIPT_DIR/$RECEPTOR_PREFIX"
OUTPUT_DIR="$BASE_DIR/vina_outputs"
RESULTS_FILE="$BASE_DIR/vina_results.csv"
ERROR_LOG="$BASE_DIR/failed_ligands.log"

mkdir -p "$PREPARED_LIG_DIR" "$OUTPUT_DIR"

# === Prepare receptor ===
echo "Preparing receptor..."
mk_prepare_receptor.py -i "$RECEPTOR" -o "$PREPARED_REC_DIR/$RECEPTOR_PREFIX" -p -v -g --box_size $SIZE --box_center $CENTER --allow_bad_res 
# mk_prepare_receptor.py --pdb "$RECEPTOR" -o "$PREPARED_REC_DIR/$RECEPTOR_PREFIX" --box_size $SIZE --box_center $CENTER --allow_bad_res

echo "Running autogrid..."
cd "$PREPARED_REC_DIR" || exit 1
autogrid4 -p "$RECEPTOR_PREFIX.gpf" -l "$RECEPTOR_PREFIX.glg"
cd -  # go back to the previous directory

# # === Initialize files ===
echo "ligand,affinity_kcal/mol" > "$RESULTS_FILE"
echo "# Failed ligands and reasons" > "$ERROR_LOG"
echo "# Format: ligand_name,step_failed,error_message" >> "$ERROR_LOG"

# # === Counters for statistics ===
total_ligands=0
successful_docking=0
failed_scrubbing=0
failed_preparation=0
failed_docking=0

# === Process each ligand ===
echo "Processing ligands..."
shopt -s nullglob
sdf_files=("$LIGAND_DIR"/*.sdf)
if [ ${#sdf_files[@]} -eq 0 ]; then
    echo "No .sdf files found in $LIGAND_DIR"
    exit 1
fi

for sdf in "${sdf_files[@]}"; do
    base=$(basename "$sdf" .sdf)
    total_ligands=$((total_ligands + 1))
    echo "Processing ligand $total_ligands: $base"

    scrubbed="$PREPARED_LIG_DIR/${base}_scrubbed.sdf"
    pdbqt="$PREPARED_LIG_DIR/${base}_scrubbed.pdbqt"
    out="$OUTPUT_DIR/${base}_out.pdbqt"

    # Step 1: Scrub ligand
    echo "  Scrubbing..."
    if ! scrub.py "$sdf" -o "$scrubbed" 2>/tmp/scrub_error_${base}.log; then
        echo "  Failed to scrub $base"
        error_msg=$(cat /tmp/scrub_error_${base}.log 2>/dev/null || echo "Unknown scrubbing error")
        echo "$base,scrubbing,$error_msg" >> "$ERROR_LOG"
        failed_scrubbing=$((failed_scrubbing + 1))
        rm -f /tmp/scrub_error_${base}.log
        continue
    fi

    # Check if scrubbed file was created and is not empty
    if [ ! -s "$scrubbed" ]; then
        echo "  Scrubbed file is empty for $base"
        echo "$base,scrubbing,Empty output file" >> "$ERROR_LOG"
        failed_scrubbing=$((failed_scrubbing + 1))
        continue
    fi

    # Step 2: Prepare ligand
    echo "  Preparing..."
    if ! mk_prepare_ligand.py -i "$scrubbed" -o "$pdbqt" 2>/tmp/prep_error_${base}.log; then
        echo "  Failed to prepare $base"
        error_msg=$(cat /tmp/prep_error_${base}.log 2>/dev/null || echo "Unknown preparation error")
        echo "$base,preparation,$error_msg" >> "$ERROR_LOG"
        failed_preparation=$((failed_preparation + 1))
        rm -f /tmp/prep_error_${base}.log
        continue
    fi

    # Check if PDBQT file was created and is not empty
    if [ ! -s "$pdbqt" ]; then
        echo "  PDBQT file is empty for $base"
        echo "$base,preparation,Empty PDBQT output" >> "$ERROR_LOG"
        failed_preparation=$((failed_preparation + 1))
        continue
    fi

    # Step 3: Run Vina
    echo "  Docking..."
    if ! vina --ligand "$pdbqt" --maps "$PREPARED_REC_DIR/$RECEPTOR_PREFIX" --scoring ad4 --exhaustiveness "$EXHAUST" --out "$out" 2>/tmp/vina_error_${base}.log; then
        echo "  Vina failed for $base"
        error_msg=$(cat /tmp/vina_error_${base}.log 2>/dev/null || echo "Unknown docking error")
        echo "$base,docking,$error_msg" >> "$ERROR_LOG"
        failed_docking=$((failed_docking + 1))
        rm -f /tmp/vina_error_${base}.log
        continue
    fi

    # if ! vina --receptor "$PREPARED_REC_DIR/$RECEPTOR_PREFIX".pdbqt --ligand "$pdbqt" \
    #    --config "$PREPARED_REC_DIR/$RECEPTOR_PREFIX".box.txt \
    #    --exhaustiveness=8 --out "$out"; then
    #     echo "Vina failed for $base"
    #     continue
    # fi

    # Step 4: Extract best binding score
    if [ -s "$out" ]; then
        score=$(grep "^REMARK VINA RESULT:" "$out" | head -1 | awk '{print $4}')
        if [ -n "$score" ]; then
            echo "$base,$score" >> "$RESULTS_FILE"
            echo "  $base: $score kcal/mol"
            successful_docking=$((successful_docking + 1))
        else
            echo "  No valid score found for $base"
            echo "$base,docking,No valid score in output" >> "$ERROR_LOG"
            failed_docking=$((failed_docking + 1))
        fi
    else
        echo "  Empty output file for $base"
        echo "$base,docking,Empty Vina output" >> "$ERROR_LOG"
        failed_docking=$((failed_docking + 1))
    fi

    # Clean up temporary error files
    rm -f /tmp/scrub_error_${base}.log /tmp/prep_error_${base}.log /tmp/vina_error_${base}.log
done

# === Print statistics ===
echo ""
echo "=== DOCKING STATISTICS ==="
echo "Total ligands processed: $total_ligands"
echo "Successful dockings: $successful_docking"
echo "Failed at scrubbing: $failed_scrubbing"
echo "Failed at preparation: $failed_preparation"
echo "Failed at docking: $failed_docking"
echo "Success rate: $(( successful_docking * 100 / total_ligands ))%"
echo ""
echo "Detailed error log saved to: $ERROR_LOG"
echo "Results saved to: $RESULTS_FILE"

# Exit with error if no ligands were successfully docked
if [ $successful_docking -eq 0 ]; then
    echo "ERROR: No ligands were successfully docked!"
    exit 1
fi

echo "All docking complete! Results saved to $RESULTS_FILE"
