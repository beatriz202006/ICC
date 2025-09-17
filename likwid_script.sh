#!/bin/bash
# Uso: ./likwid_script.sh arquivo_entrada

likwid-perfctr -C 0 -g FLOPS_DP -m ./resolveEDO < "$1" > saida_likwid.txt 2>&1

grep FP_ARITH_INST_RETIRED_SCALAR_DOUBLE saida_likwid.txt | awk -F '|' '{print "FP_ARITH_INST_RETIRED_SCALAR_DOUBLE," $4}' | sed 's/ //g'