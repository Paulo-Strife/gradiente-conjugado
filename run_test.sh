#!/bin/bash

mkdir -p resultados_perf

gcc conjugated_method.c -o conjugated_method -Wall -Wextra -O2 -lm

if [ $? -ne 0 ]; then
    echo "Erro na compilacao."
    exit 1
fi

seed=42
tol="1e-8"
max_iter=10000

eventos="cycles,instructions,cache-references,cache-misses,branches,branch-misses"

for banda in 7 27
do
    for N in 1024 4096 16384 65536 262144 1048576
    do
        echo "Rodando N=$N banda=$banda"

        sudo perf stat -e $eventos \
            ./conjugated_method $N $banda $seed $tol $max_iter \
            > resultados_perf/programa_N${N}_B${banda}.txt \
            2> resultados_perf/perf_N${N}_B${banda}.txt
    done
done

echo "Testes finalizados. Resultados salvos em resultados_perf/"