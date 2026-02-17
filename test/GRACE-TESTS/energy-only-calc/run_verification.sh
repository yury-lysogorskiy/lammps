#!/bin/bash
set -e

LMP_EXE="/home/users/lysogy36/CLionProjects/lammps-ace/build-grace/lmp"
POT_DIR="/home/users/lysogy36/CLionProjects/lammps-ace/test/GRACE-TESTS/potentials"
INPUT_TEMPLATE="/home/users/lysogy36/CLionProjects/lammps-ace/test/GRACE-TESTS/energy-only-calc/in.lammps.template"

# Array of styles to test
# Format: "STYLE_NAME POTENTIAL_FILE ELEMENTS"
STYLES=(
    "grace $POT_DIR/GRACE-2L-OMAT-medium-base-aux W Mo"
    "grace/fs $POT_DIR/FS_model.yaml W Mo"
    "grace/2layer/parallel $POT_DIR/GRACE-2L-OMAT-medium-base-aux W Mo"
    "grace/2layer/chunk $POT_DIR/GRACE-2L-OMAT-medium-base-aux W Mo"
    "grace/1layer/chunk $POT_DIR/GRACE-1L-OMAT-large-base-aux W Mo"
)

mkdir -p results

for style_info in "${STYLES[@]}"; do
    STYLE=$(echo $style_info | awk '{print $1}')
    POT=$(echo $style_info | awk '{print $2}')
    # Extract elements (everything from 3rd word onwards)
    ELEMS=$(echo $style_info | cut -d' ' -f3-)
    
    FOLDER_NAME=$(echo $STYLE | sed 's/\//_/g')
    echo "===================================================="
    echo "Testing Style: $STYLE"
    echo "===================================================="
    
    STYLE_RES_DIR="results/$FOLDER_NAME"
    mkdir -p "$STYLE_RES_DIR"

    # 1. Run Default (Energy-only optimizations active)
    echo "  Running Default..."
    RUN_DIR="$STYLE_RES_DIR/default"
    mkdir -p "$RUN_DIR"
    sed "s|POTENTIAL_STYLE|$STYLE|g; s|POTENTIAL_FILE|$POT|g; s|W Mo|$ELEMS|g" "$INPUT_TEMPLATE" > "$RUN_DIR/in.lammps"
    (
        cd "$RUN_DIR"
        $LMP_EXE -in in.lammps > screen.log 2>&1 || true
    )
    
    # 2. Run Debug (Energy-only optimizations disabled)
    echo "  Running Debug (no optimizations)..."
    RUN_DIR="$STYLE_RES_DIR/debug"
    mkdir -p "$RUN_DIR"
    sed "s|POTENTIAL_STYLE|$STYLE debug_no_energy_only_calc|g; s|POTENTIAL_FILE|$POT|g; s|W Mo|$ELEMS|g" "$INPUT_TEMPLATE" > "$RUN_DIR/in.lammps"
    (
        cd "$RUN_DIR"
        $LMP_EXE -in in.lammps > screen.log 2>&1 || true
    )
    
    # Extraction and Verification
    LOG_DEFAULT="$STYLE_RES_DIR/default/log.lammps"
    LOG_DEBUG="$STYLE_RES_DIR/debug/log.lammps"

    if [[ -f "$LOG_DEFAULT" && -f "$LOG_DEBUG" ]]; then
        grep -A 10 "Step" "$LOG_DEFAULT" | grep -v "Step" | head -n 11 > "$STYLE_RES_DIR/thermo.default"
        grep -A 10 "Step" "$LOG_DEBUG" | grep -v "Step" | head -n 11 > "$STYLE_RES_DIR/thermo.debug"
        
        echo "Verification for $STYLE:"
        if diff "$STYLE_RES_DIR/thermo.default" "$STYLE_RES_DIR/thermo.debug" > /dev/null; then
            echo "    [SUCCESS] Thermodynamic data matches."
        else
            echo "    [FAILED] Thermodynamic data MISMATCH!"
            diff -u "$STYLE_RES_DIR/thermo.default" "$STYLE_RES_DIR/thermo.debug"
        fi
        
        T_DEFAULT=$(grep "Loop time of" "$STYLE_RES_DIR/default/screen.log" | awk '{print $4}' || echo "N/A")
        T_DEBUG=$(grep "Loop time of" "$STYLE_RES_DIR/debug/screen.log" | awk '{print $4}' || echo "N/A")
        echo "    Execution Times: Default: ${T_DEFAULT}s, Debug: ${T_DEBUG}s"
    else
        echo "    [ERROR] Logs missing for $STYLE. Check screen.log in results/$FOLDER_NAME/"
    fi
done

echo "===================================================="
echo "Verification Complete."
