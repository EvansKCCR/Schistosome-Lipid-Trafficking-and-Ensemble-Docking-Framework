#!/bin/bash

RECEPTOR1="receptor_grid/receptor_rigid.pdbqt"
RECEPTOR2="receptor_grid/receptor_flex.pdbqt"
BOX="receptor_grid/receptor_box.txt"

OUTDIR="vina_results"
LOGDIR="vina_logs"

mkdir -p $OUTDIR $LOGDIR

# Extract box parameters
CENTER_X=$(grep center_x $BOX | awk '{print $3}')
CENTER_Y=$(grep center_y $BOX | awk '{print $3}')
CENTER_Z=$(grep center_z $BOX | awk '{print $3}')

SIZE_X=$(grep size_x $BOX | awk '{print $3}')
SIZE_Y=$(grep size_y $BOX | awk '{print $3}')
SIZE_Z=$(grep size_z $BOX | awk '{print $3}')

# ==============================
# DOCKING LOOP
# ==============================

for lig in ligand_pdbqt/*.pdbqt; do

    base=$(basename "$lig" .pdbqt)

    echo "🔄 Running Vina for $base"

    LIG_OUTDIR=${OUTDIR}/${base}
    mkdir -p $LIG_OUTDIR

    vina \
        --receptor $RECEPTOR1 \
        --flex $RECEPTOR2 \
        --ligand "$lig" \
        --center_x $CENTER_X \
        --center_y $CENTER_Y \
        --center_z $CENTER_Z \
        --size_x $SIZE_X \
        --size_y $SIZE_Y \
        --size_z $SIZE_Z \
        --exhaustiveness 64 \
        --num_modes 50 \
        --energy_range 3 \
        --seed 42 \
        --out ${LIG_OUTDIR}/${base}_flex_out.pdbqt \
        > ${LOGDIR}/${base}.log 2>&1

    echo "✅ Completed docking: $base"

done

echo "🎉 All Vina refinements completed"