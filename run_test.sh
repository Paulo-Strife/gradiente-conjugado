#!/bin/bash

gcc conjugated_method.c -o conjugated_method -Wall -Wextra -O2 -lm

if [ $? -ne 0 ]; then
    echo "Erro na compilacao. O executavel nao foi gerado."
    exit 1
fi

seed=42
tol="1e-8"
max_iter=10000

for banda in 7 27
do
    for N in 1024 4096 16384 65536 262144 1048576
    do
        echo "========================================"
        echo "Rodando N=$N banda=$banda"
        echo "========================================"

        ./conjugated_method $N $banda $seed $tol $max_iter

        echo ""
    done
done