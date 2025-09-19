# Beatriz Pontes Camargo
# GRR 20242966
./resolveEDO < "$1"
#!/bin/bash
# Uso: ./likwid_script.sh arquivo_entrada  => exemplo: ./likwid_script.sh teste.dat
# ./resolveEDO < "$1" → executa o programa usando como entrada o arquivo recebido no argumento ($1)
# Saída padrão e de erro vão para o arquivo saida_likwid.txt
likwid-perfctr -C 0 -g FLOPS_DP -m ./resolveEDO < "$1" > saida_likwid.txt 2>&1


#Filtra o relatório do LIKWID para pegar apenas a linha com o contador:
# awk -F '|' '{print "FP_ARITH_INST_RETIRED_SCALAR_DOUBLE," $4}'
#   -F '|' → define o separador de campos como "|"
#   $4     → pega o quarto campo da linha
# sed 's/ //g' remove todos os espaços em branco da saída final
grep FP_ARITH_INST_RETIRED_SCALAR_DOUBLE saida_likwid.txt | awk -F '|' '{print "FP_ARITH_INST_RETIRED_SCALAR_DOUBLE," $4}' | sed 's/ //g'
