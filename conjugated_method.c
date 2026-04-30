#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct no {
    float valor;
    int coluna;
    struct no* prox;
} NO;

typedef NO* PONT;

typedef struct {
    PONT* A;
    int linhas;
    int colunas;
    int banda;
    int semiBanda;
} MATRIZ;

bool dentroBanda(MATRIZ* m, int lin, int col) {
    return abs(lin-col) <= m->semiBanda; 
}

bool inicializaMatriz(MATRIZ* m, int lin, int col, int banda) {
    int i;

    if (banda != 7 && banda != 27)
        return false;
    
    m->linhas = lin;
    m->colunas = col;
    m->banda = banda;
    m->semiBanda = (banda - 1) / 2;

    m->A = (PONT *) malloc(lin* sizeof(PONT));
    if (m->A == NULL)
        return false;

    for (i = 0; i<lin; i++)
        m->A[i] = NULL;

    return true;
}

bool atribuiMatriz(MATRIZ* m, int lin, int col, float val) {
    if (lin < 0 || lin >= m->linhas || col < 0 || col >= m->colunas) {
        return false;
    }

    if (!dentroBanda(m, lin, col)) {
        if (val != 0.0f) {
            printf("Posicao (%d, %d) fora da banda permitida da matriz de %d bandas.\n",
                   lin, col, m->banda);
            return false;
        } else {
            return true;
        }
    }

    PONT ant = NULL;
    PONT atual = m->A[lin];

    while (atual != NULL && atual->coluna < col) {
        ant = atual;
        atual = atual->prox;
    }

    if (atual != NULL && atual->coluna == col) {
        if (val == 0.0f) {
            if (ant == NULL) {
                m->A[lin] = atual->prox;
            } else {
                ant->prox = atual->prox;
            }
            free(atual);
        } else {
            atual->valor = val;
        }
        return true;
    }

    if (val != 0.0f) {
        PONT novo = (PONT) malloc(sizeof(NO));
        if (novo == NULL) {
            return false;
        }

        novo->coluna = col;
        novo->valor = val;
        novo->prox = atual;

        if (ant == NULL) {
            m->A[lin] = novo;
        } else {
            ant->prox = novo;
        }
    }

    return true;
}

float valorMatriz(MATRIZ* m, int lin, int col) {
    if (lin < 0 || lin >= m->linhas || col < 0 || col >= m->colunas) {
        return 0.0f;
    }

    if (!dentroBanda(m, lin, col)) {
        return 0.0f;
    }

    PONT atual = m->A[lin];

    while (atual != NULL && atual->coluna < col) {
        atual = atual->prox;
    }

    if (atual != NULL && atual->coluna == col) {
        return atual->valor;
    }

    return 0.0f;
}

void imprimirMatrizEsparsa(MATRIZ* m) {
    int i;
    for (i = 0; i < m->linhas; i++) {
        PONT atual = m->A[i];
        printf("Linha %d: ", i);
        while (atual != NULL) {
            printf("(%d, %.1f) ", atual->coluna, atual->valor);
            atual = atual->prox;
        }
        printf("\n");
    }
}

void apagaMatriz(MATRIZ* m) {
    int i;
    for (i = 0; i<m->linhas; i++) {
        PONT atual = m->A[i];
        while (atual != NULL) {
            PONT temp = atual;
            atual = atual->prox;
            free(temp);
        }
    }
    free(m->A);
    m->A = NULL;
}

void preencheBanda(MATRIZ* m) {
    int i, desloc;

    for (i = 0; i<m->linhas; i++) {
        atribuiMatriz(m, i, i, 10.0f);

        for (desloc = 1; desloc <= m->semiBanda; desloc++) {
            if (i + desloc < m->colunas) {
                atribuiMatriz(m, i, i + desloc, 1.0f * desloc);
            }
            if (i - desloc >= 0) {
                atribuiMatriz(m, i, i - desloc, 1.0f * desloc);
            }
        }
    }
}

int main() {
    MATRIZ m;

    if (!inicializaMatriz(&m, 10, 10, 7)) {
        printf("Erro ao incializar o programa");
        return 1;
    }
    printf("bandas inicializadas %d. \n", m.banda);
    printf("Semi-Banda = %d\n\n", m.semiBanda);

    preencheBanda(&m);
    imprimirMatrizEsparsa(&m);
    return 0;
}