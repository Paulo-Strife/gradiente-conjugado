#!/bin/bash

mkdir -p resultados

SRC="conjugated_method.c"
BIN="conjugated_method"

gcc $SRC -o $BIN -Wall -Wextra -O0 -lm

if [ $? -ne 0 ]; then
    echo "Erro na compilacao."
    exit 1
fi

seed=42
tol="1e-8"
max_iter=10000
repeticoes=5

eventos="cycles,instructions,cache-references,cache-misses,branches,branch-misses"

echo "N,banda,rep,tempo_s,memoria_kb,cycles,instructions,cache_references,cache_misses,branches,branch_misses" > resultados/benchmark_atual.csv

for banda in 7 27
do
    for N in 1024 4096 16384 65536 262144 1048576
    do
        for rep in $(seq 1 $repeticoes)
        do
            echo "Rodando N=$N banda=$banda rep=$rep"

            sudo perf stat -x, -e $eventos \
                /usr/bin/time -v ./$BIN $N $banda $seed $tol $max_iter \
                > resultados/temp_programa.txt \
                2> resultados/temp_metricas.txt

            tempo_s=$(grep "Tempo total CG" resultados/temp_programa.txt | awk -F ':' '{print $2}' | sed 's/s//g' | xargs)

            memoria_kb=$(grep "Maximum resident set size" resultados/temp_metricas.txt | awk -F ':' '{print $2}' | xargs)

            cycles=$(grep ",cycles," resultados/temp_metricas.txt | head -n 1 | awk -F ',' '{print $1}' | tr -d ' ')
            instructions=$(grep ",instructions," resultados/temp_metricas.txt | head -n 1 | awk -F ',' '{print $1}' | tr -d ' ')
            cache_references=$(grep ",cache-references," resultados/temp_metricas.txt | head -n 1 | awk -F ',' '{print $1}' | tr -d ' ')
            cache_misses=$(grep ",cache-misses," resultados/temp_metricas.txt | head -n 1 | awk -F ',' '{print $1}' | tr -d ' ')
            branches=$(grep ",branches," resultados/temp_metricas.txt | head -n 1 | awk -F ',' '{print $1}' | tr -d ' ')
            branch_misses=$(grep ",branch-misses," resultados/temp_metricas.txt | head -n 1 | awk -F ',' '{print $1}' | tr -d ' ')

            echo "$N,$banda,$rep,$tempo_s,$memoria_kb,$cycles,$instructions,$cache_references,$cache_misses,$branches,$branch_misses" >> resultados/benchmark_atual.csv
        done
    done
done

rm -f resultados/temp_programa.txt resultados/temp_metricas.txt

echo "Arquivo gerado: resultados/benchmark_atual.csv"