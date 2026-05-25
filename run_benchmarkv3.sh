#!/bin/bash

mkdir -p resultados

SRC_BASELINE="conjugated_method.c"
SRC_OPT="conjugated_methodv3.c"

BIN_BASELINE="baseline"
BIN_OPT="opt_unroll"

FLAGS="-Wall -Wextra -O0"

gcc $SRC_BASELINE -o $BIN_BASELINE $FLAGS -lm

if [ $? -ne 0 ]; then
    echo "Erro ao compilar baseline."
    exit 1
fi

gcc $SRC_OPT -o $BIN_OPT $FLAGS -lm

if [ $? -ne 0 ]; then
    echo "Erro ao compilar opt_unroll."
    exit 1
fi

seed=42
tol="1e-8"
max_iter=10000
repeticoes=5

eventos="cycles,instructions,cache-references,cache-misses,branches,branch-misses"

echo "versao,N,banda,rep,tempo_s,memoria_kb,cycles,instructions,cache_references,cache_misses,branches,branch_misses" > resultados/benchmark_unroll.csv

pegar_evento() {
    arquivo="$1"
    evento="$2"

    awk -F',' -v evento="$evento" '
    {
        valor=$1
        nome=$3

        gsub(/^[ \t]+|[ \t]+$/, "", valor)
        gsub(/^[ \t]+|[ \t]+$/, "", nome)

        if (valor == "<not counted>") {
            next
        }

        if (nome == evento || nome ~ "/" evento "/") {
            gsub(/[^0-9.]/, "", valor)
            print valor
            exit
        }
    }
    ' "$arquivo"
}

for versao in baseline opt_unroll
do
    if [ "$versao" = "baseline" ]; then
        bin="./$BIN_BASELINE"
    else
        bin="./$BIN_OPT"
    fi

    for banda in 7 27
    do
        for N in 1024 4096 16384 65536 262144 1048576
        do
            for rep in $(seq 1 $repeticoes)
            do
                echo "Rodando versao=$versao N=$N banda=$banda rep=$rep"

                sudo perf stat -x, -o resultados/temp_perf.csv -e $eventos -- \
                    /usr/bin/time -v -o resultados/temp_time.txt \
                    $bin $N $banda $seed $tol $max_iter \
                    > resultados/temp_programa.txt

                tempo_s=$(grep "Tempo total CG" resultados/temp_programa.txt | awk -F ':' '{print $2}' | sed 's/s//g' | xargs)
                memoria_kb=$(grep "Maximum resident set size" resultados/temp_time.txt | awk -F ':' '{print $2}' | xargs)

                cycles=$(pegar_evento resultados/temp_perf.csv "cycles")
                instructions=$(pegar_evento resultados/temp_perf.csv "instructions")
                cache_references=$(pegar_evento resultados/temp_perf.csv "cache-references")
                cache_misses=$(pegar_evento resultados/temp_perf.csv "cache-misses")
                branches=$(pegar_evento resultados/temp_perf.csv "branches")
                branch_misses=$(pegar_evento resultados/temp_perf.csv "branch-misses")

                echo "$versao,$N,$banda,$rep,$tempo_s,$memoria_kb,$cycles,$instructions,$cache_references,$cache_misses,$branches,$branch_misses" >> resultados/benchmark_unroll.csv
            done
        done
    done
done

rm -f resultados/temp_programa.txt resultados/temp_perf.csv resultados/temp_time.txt

echo "Arquivo gerado: resultados/benchmark_unroll.csv"