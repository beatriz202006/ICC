#!/bin/bash

# Programas
v1=./v1/cgSolver
v2=./v2/cgSolver
likwid=/home/soft/likwid/bin/likwid-perfctr

# Tamanhos
Ns=(32 64 128 256)

# Grupos válidos do LIKWID
GRUPOS=("FLOPS_DP" "L2" "L3" "DATA")

echo "Starting tests: $(date)"
echo "v1 = $v1"
echo "v2 = $v2"
echo "likwid = $likwid"

for N in "${Ns[@]}"; do
    echo ""
    echo "-------------------------------------------"
    echo " Testando N = $N"
    echo "-------------------------------------------"

    for grp in "${GRUPOS[@]}"; do
        echo "-- Executando v1  (grupo $grp) --"

        $likwid -C 0 -g "$grp" -m $v1 <<< "$N 7 0.0 25 1e-6" \
            > saida_v1_${N}_${grp}.txt

        echo "-- Executando v2  (grupo $grp) --"

        $likwid -C 0 -g "$grp" -m $v2 <<< "$N 7 0.0 25 1e-6" \
            > saida_v2_${N}_${grp}.txt
    done
done

