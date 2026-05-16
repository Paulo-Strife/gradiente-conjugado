#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Colocando as variáveis principais do código
typedef struct {
    int n;                          // tamanho da matriz n x n
    int banda;                      // numero de bandas
    int meia_banda;                 // quantidade de bandas acima da diagonal principal
    double **diag;                  // contador da diagonal principal (não sei o que estou fazendo kkkkk)
} bandaMatriz;

typedef struct {
    double time_matvec;
    double time_ponto;
    double time_att;
    int iteracoes;
    double final_residual;
} GCstatus;

/*
A função temp_calc é uma função que usa um dos relógios do sistema para medir 
o tempo de execução de um código em segundos.
*/
double temp_calc() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

/*Essa função é necessária para gerar os números aleatórios entre o valor mínimo e o valor máximo
Usaremos o valor criado dentro desta função para popular as matrizes com o valor aleatórios necessário */
double rand_double(double min, double max) {
    return min + (max - min) * ((double) rand() / (double) RAND_MAX);
}

bandaMatriz *cria_matriz(int n, int banda) {
    bandaMatriz *A = (bandaMatriz *) malloc(sizeof(bandaMatriz));

    if (A == NULL) {
        printf("Erro ao formar a matriz.\n");
        exit(1);
    }

    A->n = n;
    A->banda = banda;
    A->meia_banda = (banda - 1) / 2;

    A->diag = (double **) malloc((A->meia_banda + 1) * sizeof(double *));

    if (A->diag == NULL) {
        printf("Erro ao formar a diagonal e alocar num vetor\n");
        exit(1);
    }

    for (int d = 0; d <= A->meia_banda; d++) {
        A->diag[d] = (double *) calloc(n, sizeof(double));

        if (A->diag[d] == NULL) {
            printf("Erro ao alocar a diagonal %d.\n", d);
            exit(1);
        }
    }
    return A;
}

// Libera a memória da matriz

void libera_memoria(bandaMatriz *A) {
    if (A == NULL) return;

    for (int d = 0; d <= A->meia_banda; d++) {
        free(A->diag[d]);
    }

    free(A->diag);
    free(A);
}

/*Aqui vamos gerar uma matriz simétrica positiva
    - Aqui geramos valores aleatórios pequenos nas bandas superiores
    - Depois somamos a diagonal principal.
*/

void gerar_banda_espalhada_matriz(bandaMatriz *A, unsigned int seed) {
    srand(seed);

    int n = A->n;
    int mb = A->meia_banda;

    double *seq_soma = (double *) calloc(n, sizeof(double));

    if (seq_soma == NULL) {
        printf("erro ao colocar o seq_sum\n");
        exit(1);
    }

    /*
    Gera as diagonias superiores
    */
   for (int d = 1; d <= mb; d++) {
    for (int i = 0; i < n - d; i++) {
        double valor = rand_double(-0.5, 0.5);

        A->diag[d][i] =  valor;

        seq_soma[i] += fabs(valor);
        seq_soma[i + d] += fabs(valor);
    }
   }

   /*
   Diagonal Principal, fiz ela ser maior que a soma dos outros elementos da linha
   */
  for (int i = 0; i < n; i++) {
    A->diag[0][i] = seq_soma[i] + rand_double(1.0, 2.0);
  }

  free(seq_soma);
}

// essa função é responsável por multiplicar as matrizes como no método do gradiente conjugado - básicamente
// um A * X = Y.
void multiplicacao_matrizes(const bandaMatriz *A, const double *x, double *y) {
    int n = A->n;
    int mb = A->meia_banda;

    for (int i = 0; i < n; i++) {
        y[i] = A->diag[0][i] * x[i];
    }

    for (int d = 1; d<=mb; d++) {
        for (int i = 0; i < n - d; i++) {
            double valor = A->diag[d][i];
            
            // parte superior
            y[i] += valor * x[i + d];
            // parte inferior
            y[i + d] += valor * x[i];
        }
    }
}

// produto de a*T = b
double produto(const double *a, const double *b, int n) {
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }

    return sum;
}

double vetor_norma(const double *x, int n) {
    return sqrt(produto(x, x, n));
}

// função para preencher o vetor com o valor constante
void preenche_vetor(double *v, int n, double valor) {
    for (int i = 0; i < n; i++) {
        v[i] = valor;
    }
}

/*
Agora começaremos a aplicar o método do gradiente conjugado

A x = b
*/
void print_vetor(const char *nome, const double *v, int n) {
    printf("%s = [", nome);

    for (int i = 0; i < n; i++) {
        printf("%.8f", v[i]);

        if (i < n - 1) {
            printf(", ");
        }
    }

    printf("]\n");
}

double pega_valor_matriz(const bandaMatriz *A, int i, int j) {
    int diff = j - i;

    if (diff < 0) {
        diff = -diff;
    }

    if (diff > A->meia_banda) {
        return 0.0;
    }

    if (j >= i) {
        return A->diag[diff][i];
    } else {
        return A->diag[diff][j];
    }
}

void print_matriz_completa(const bandaMatriz *A) {
    for (int i = 0; i < A->n; i++){
        printf("[ ");

        for (int j = 0; j < A->n; j++) {
            printf("%.12f ", pega_valor_matriz(A, i, j));
        }

        printf("]\n");
    }

    printf("\n");
}

int gradiente_conjugado(
    const bandaMatriz *A,
    const double *b,
    double *x,
    int iter,
    double tol,
    GCstatus *stats,
    int debug
) {
    int n = A->n;

    if (debug) {
        printf("\n====================== MATRIZ DO SISTEMA ======================\n");
        print_matriz_completa(A);
    }

    double *r = (double *) malloc(n * sizeof(double));
    double *p = (double *) malloc(n * sizeof(double));
    double *Ap = (double *) malloc(n * sizeof(double));

    if (r == NULL || p == NULL || Ap == NULL) {
        printf("Erro ao alocar vetores");
        exit(1);
    }

    stats->time_matvec = 0.0;
    stats->time_ponto = 0.0;
    stats->time_att = 0.0;
    stats->iteracoes = 0;
    stats->final_residual = 0.0;

    double t0 = temp_calc();
    multiplicacao_matrizes(A, x, Ap);
    double t1 = temp_calc();

    stats->time_matvec += t1 - t0;

    t0 = temp_calc();
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - Ap[i];
        p[i] = r[i];
    }
    t1 = temp_calc();

    stats->time_att += t1 - t0;

    if (debug) {
        printf("\n========== INICIO DO GRADIENTE CONJUGADO ==========\n");
        print_vetor("x inicial", x, n);
        print_vetor("b", b, n);
        print_vetor("Ax inicial", Ap, n);
        print_vetor("r inicial = b - Ax", r, n);
        print_vetor("p inicial = r", p, n);
    }

    t0 = temp_calc();
    double rs_antigo = produto(r, r, n);
    t1 = temp_calc();

    stats->time_ponto += t1 - t0;

    double residuo_inicial = sqrt(rs_antigo);

    if (debug) {
        printf("\nrs_antigo = r^T r = %.12f\n", rs_antigo);
        printf("residuo inicial = sqrt(rs_antigo) = %.12f\n", residuo_inicial);
    }

    if (residuo_inicial < tol) {
        stats->final_residual = residuo_inicial;
        free(r);
        free(p);
        free(Ap);
        return 0;
    }

    for (int k = 0; k < iter; k++) {
        if (debug) {
            printf("\n========================================\n");
            printf("ITERACAO %d\n", k + 1);
            printf("========================================\n");
            print_vetor("x atual", x, n);
            print_vetor("r atual", r, n);
            print_vetor("p atual", p, n);
            printf("rs_antigo = %.12f\n", rs_antigo);
        }

        t0 = temp_calc();
        multiplicacao_matrizes(A, p, Ap);
        t1 = temp_calc();

        stats->time_matvec += t1 - t0;

        if (debug) {
            printf("\n1) Ap = A * p\n");
            print_vetor("Ap", Ap, n);
        }

        t0 = temp_calc();
        double pAp = produto(p, Ap, n);
        t1 = temp_calc();

        stats->time_ponto += t1 - t0;

        if (debug) {
            printf("\n2) pAp = p^T Ap\n");
            printf("pAp = %.12f\n", pAp);
        }

        if (fabs(pAp) < 1e-30) {
            printf("Erro numerico, p^T A p é muito pequeno\n");
            break;
        }

        double alpha = rs_antigo / pAp;

        if (debug) {
            printf("\n3) alpha = rs_antigo / pAp\n");
            printf("alpha = %.12f / %.12f = %.12f\n", rs_antigo, pAp, alpha);
        }

        if (debug) {
            printf("\n4) Atualizando x e r\n");

            for (int i = 0; i < n; i++) {
                double novo_x = x[i] + alpha * p[i];
                double novo_r = r[i] - alpha * Ap[i];

                printf("x[%d] = %.12f + %.12f * %.12f = %.12f\n",
                       i, x[i], alpha, p[i], novo_x);
                printf("r[%d] = %.12f - %.12f * %.12f = %.12f\n",
                       i, r[i], alpha, Ap[i], novo_r);
            }
        }

        t0 = temp_calc();
        for (int i = 0; i < n; i++) {
            x[i] = x[i] + alpha * p[i];
            r[i] = r[i] - alpha * Ap[i];
        }
        t1 = temp_calc();

        stats->time_att += t1 - t0;

        if (debug) {
            print_vetor("x novo", x, n);
            print_vetor("r novo", r, n);
        }

        t0 = temp_calc();
        double rs_novo = produto(r, r, n);
        t1 = temp_calc();

        stats->time_ponto += t1 - t0;

        double residual = sqrt(rs_novo);

        stats->iteracoes = k + 1;
        stats->final_residual = residual;

        if (debug) {
            printf("\n5) Novo residuo\n");
            printf("rs_novo = r^T r = %.12f\n", rs_novo);
            printf("residual = sqrt(rs_novo) = %.12f\n", residual);
        }

        if (residual < tol) {
            if (debug) {
                printf("\nCriterio de parada atingido\n");
                printf("residual %.12f < tol %.12f\n", residual, tol);
            }
            break;
        }

        double beta = rs_novo / rs_antigo;

        if (debug) {
            printf("\n6) beta = rs_novo / rs_antigo\n");
            printf("beta = %.12f / %.12f = %.12f\n", rs_novo, rs_antigo, beta);
        }

        if (debug) {
            printf("\n7) Atualizando p\n");

            for (int i = 0; i < n; i++) {
                double novo_p = r[i] + beta * p[i];

                printf("p[%d] = %.12f + %.12f * %.12f = %.12f\n",
                       i, r[i], beta, p[i], novo_p);
            }
        }

        t0 = temp_calc();
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }
        t1 = temp_calc();

        stats->time_att += t1 - t0;

        if (debug) {
            print_vetor("p novo", p, n);
        }

        rs_antigo = rs_novo;

        if (debug) {
            printf("\nFim da iteracao %d\n", k + 1);
            printf("rs_antigo recebe rs_novo: %.12f\n", rs_antigo);
        }
    }

    free(r);
    free(p);
    free(Ap);

    return stats->iteracoes;
}

double erro_max(const double *x, const double *x_true, int n) {
    double max_err = 0.0;

    for (int i = 0; i < n; i++) {
        double err = fabs(x[i] - x_true[i]);

        if (err > max_err) {
            max_err = err;
        }
    }

    return max_err;
}

void free_matrix(bandaMatriz *A) {
    if (A == NULL) {
        return;
    }

    for (int d = 0; d<= A->meia_banda; d++) {
        free(A->diag[d]);
    }

    free(A->diag);
    free(A);
}

int main (int args, char **argv) {
    if (args < 6) {
        printf("Uso:\n");
        printf("%s <n> <bandas> <seed> <tol> <max_iter>\n\n", argv[0]);
        printf("Exemplos:\n");
        printf("%s 1024 7 42 1e-8", argv[0]);
        printf("%s 4096 27 42 1e-8", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int banda = atoi(argv[2]);
    unsigned int seed = (unsigned int) atoi(argv[3]);
    double tol = atof(argv[4]);
    int max_iter = atoi(argv[5]);

    if (banda != 7 && banda != 27) {
        printf("Erro, banda deve ser 7 ou 27");
        return 1;
    }

    printf("========================================\n");
    printf("Gradiente Conjugado - Matriz por bandas\n");
    printf("n         = %d\n", n);
    printf("bands     = %d\n", banda);
    printf("seed      = %u\n", seed);
    printf("tol       = %.2e\n", tol);
    printf("max_iter  = %d\n", max_iter);
    printf("========================================\n");

    double total_comeco = temp_calc();

    // cria e gera a matriz

    double t0 = temp_calc();

    bandaMatriz *A = cria_matriz(n, banda);
    gerar_banda_espalhada_matriz(A, seed);

    double t1 = temp_calc();

    double time_generation = t1 - t0;

    double *x_true = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *x =  (double *) malloc(n * sizeof(double));

    if (x_true == NULL || b == NULL || x == NULL) {
        printf("Erro ao alocar vetores");
        exit(1);
    }

    preenche_vetor(x_true, n, 1.0);
    preenche_vetor(x, n, 0.0);

    t0 = temp_calc();
    multiplicacao_matrizes(A, x_true, b);
    t1 = temp_calc();

    double time_generate_b = t1 - t0;

    GCstatus status;

    t0 =  temp_calc();
    int debug = (n <= 5);
    gradiente_conjugado(A, b, x, max_iter, tol, &status, debug);
    t1 = temp_calc();

    double tempo_gc = t1 - t0;

    double err = erro_max(x, x_true, n);

    double total_final = temp_calc();

    printf("\nResultado:\n");
    printf("Iteracoes                  : %d\n", status.iteracoes);
    printf("Residuo final              : %.12e\n", status.final_residual);
    printf("Erro maximo contra x_true  : %.12e\n", err);

    printf("\nTempos:\n");
    printf("Geracao da matriz          : %.6f s\n", time_generation);
    printf("Geracao do vetor b         : %.6f s\n", time_generate_b);
    printf("Tempo total CG             : %.6f s\n", tempo_gc);
    printf("  Tempo em matvec          : %.6f s\n", status.time_matvec);
    printf("  Tempo em dot products    : %.6f s\n", status.time_ponto);
    printf("  Tempo em updates vetores : %.6f s\n", status.time_att);
    printf("Tempo total do programa    : %.6f s\n", total_final - total_comeco);

    // liberar memória

    free(x_true);
    free(b);
    free(x);
    free_matrix(A);

    return 0;
}