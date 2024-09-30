/*
Resumo do Código
Definições e Inicializações: O código começa definindo constantes e inicializando variáveis para simulação de dinâmica de fluidos 
usando o método Lattice Boltzmann.

Cálculo de Equilíbrio e Colisão: Dentro do laço de tempo, ele calcula a distribuição de equilíbrio e atualiza as distribuições de partículas.

Streaming: Realiza o movimento das partículas (streaming) nas direções definidas.

Condições de Contorno: Define como as partículas interagem nas bordas do domínio, garantindo que as condições desejadas sejam respeitadas.

Cálculo da Densidade e Plotagem: Por fim, calcula a densidade e prepara dados para visualização ou análise.
*/

#include <stdio.h>    // biblioteca padrão de entrada e saída
#include <stdlib.h>   // biblioteca para alocação de memória
#include <string.h>   // biblioteca para manipulação de strings

// constantes para o tamanho das matrizes e número de direções, nesse caso D2Q9
#define M 40         // Número de pontos em x
#define N 40         // Número de pontos em y
#define Q 9          // Número de direções

// função auxiliar para circular a propagação das matrizes, substituindo o circshift do MATLAB
void streaming_direction(double f[M][N][Q], double temp[M][N][Q], int dx, int dy, int k) {
    // dx e dy representam o deslocamento em x e y respectivamente, que indicam de onde a informação deve ser trazida/propagada para o ponto (i,j)
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            int new_i = i + dx;
            int new_j = j + dy;

            // verifica se os índices estão dentro dos limites da matriz
            if (new_i >= 0 && new_i < M && new_j >= 0 && new_j < N) {
                f[i][j][k] = temp[new_i][new_j][k];
            }
        }
    }
}


// Função principal
int main() {
    // Definindo variáveis
    double xl = 1.0, yl = 1.0;    // Dimensões do domínio
    double dx = xl / (M - 1.0);    // Passo em x
    double dy = yl / (N - 1.0);    // Passo em y
    double w0 = 4.0 / 9.0;         // Peso para a direção de repouso, o valor 4.0 / 9.0 é o peso associado à direção estacionária (sem movimento) no esquema de discretização D2Q9 A.A. Mohamad pg44
    // Matrizes para distribuição
    double f0[M][N],  // função de distribuição na direção estacionária (direção sem movimento)
    f0eq[M][N],       // função de distribuição de equilíbrio.
    f[M][N][Q],       // funções de distribuição (para cada uma das 9 direções do modelo D2Q9)
    feq[M][N][Q];     // funções de distribuição de equilíbrio.
    
    double rho[M][N], x[M], y[N], fluxq[M], flux[M]; // Matrizes e vetores auxiliares
    double Tm[M], Z[N][M], w[Q];   // Matrizes e vetores para resultados e pesos
    double alpha = 0.75;            // Parâmetro para o modelo (0,25)
    double omega = 1.0 / (3.0 * alpha + 0.5); // Relaxação
    double twall = 1.0;             // Condição da parede
    int nstep = 100;                // Número de passos de tempo (400)

    // Inicializando arrays com zeros
    memset(f0, 0, sizeof(f0));     
    memset(f, 0, sizeof(f));       
    memset(feq, 0, sizeof(feq));   
    memset(f0eq, 0, sizeof(f0eq)); 
    memset(rho, 0, sizeof(rho));   
    memset(fluxq, 0, sizeof(fluxq));
    memset(flux, 0, sizeof(flux));  
    memset(Tm, 0, sizeof(Tm));      
    memset(Z, 0, sizeof(Z));        
    memset(w, 0, sizeof(w));        
    {/*função MEMSET
    Resumo dos Parâmetros:
    1 (ptr): ponteiro para o bloco de memória que você quer preencher.
    2 (value): valor a ser usado para preencher cada byte do bloco de memória.
    3 (num): número de bytes a serem preenchidos com o valor especificado.*/
    }
    
    //Preenchendo valores das coordenadas 
    for (int i = 0; i < M; i++) {
        x[i] = i * dx;
    }
    
    for (int j = 0; j < N; j++) {
        y[j] = j * dy;
    }
    
    // Inicializando os pesos
    for (int k = 0; k < 4; k++) {
        w[k] = 1.0 / 9.0;      // Pesos para as quatro direções cardeais, w(1) = w(2) = w(3) = w(4) = 1/9; (muhamad)
    }
    for (int k = 4; k < Q; k++) {
        w[k] = 1.0 / 36.0;     // Pesos para as outras direções diagonais, w(5) = w(6) = w(7) = w(8) = 1/36; (muhamad)
    }

    // Colisão   fk(x, y,t + deltat) = fk(x,y,t)[1−delta]+ omega*feqk (x,y,t)
    for (int kk = 0; kk < nstep; kk++) { // Laço de tempo
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < M; i++) {
                // atualizações na direção estacionaria
                f0eq[i][j] = w0 * rho[i][j]; // função de distribuição de equilíbrio local da densidade 
                f0[i][j] = (1.0 - omega) * f0[i][j] + omega * f0eq[i][j]; // função de distribuição com base na função de equilíbrio f0eq[i][j] e no parâmetro de relaxação omega (BGK).
                
                // atualização nas outras direções
                for (int k = 0; k < Q; k++) {
                    feq[i][j][k] = w[k] * rho[i][j]; // função de distribuição de equilíbrio local da densidade 
                    f[i][j][k] = (1.0 - omega) * f[i][j][k] + omega * feq[i][j][k]; // Atualiza a função de distribuição f0[i][j] com base na função de equilíbrio f0eq[i][j] e no parâmetro de relaxação omega (BGK).
                }
            }
        }
    }

    // propagação
    double temp[M][N][Q]; // Array temporário para armazenar os dados
    memcpy(temp, f, sizeof(f)); // Copia os dados de f para temp

    //propagação para cada uma das direções 1 a 8 (colocar dentro de um laço for)
    streaming_direction(f, temp, +1,  0, 0); // Direção 1 (dx=+1, dy=0)
    streaming_direction(f, temp,  0, -1, 1); // Direção 2 (dx=0, dy=-1)
    streaming_direction(f, temp, -1,  0, 2); // Direção 3 (dx=-1, dy=0)
    streaming_direction(f, temp,  0, +1, 3); // Direção 4 (dx=0, dy=+1)
    streaming_direction(f, temp, +1, +1, 4); // Direção 5 (dx=+1, dy=+1)
    streaming_direction(f, temp, -1, +1, 5); // Direção 6 (dx=-1, dy=+1)
    streaming_direction(f, temp, -1, -1, 6); // Direção 7 (dx=-1, dy=-1)
    streaming_direction(f, temp, +1, -1, 7); // Direção 8 (dx=+1, dy=-1)

    // Condições de contorno
    // Condição de contorno esquerda, twall=1.0
    for (int j = 0; j < N; j++) {
        f[0][j][0] = w[0] * twall + w[2] * twall - f[0][j][2]; // Ajusta direção 1
        f[0][j][4] = w[4] * twall + w[6] * twall - f[0][j][6]; // Ajusta direção 5
        f[0][j][7] = w[7] * twall + w[5] * twall - f[0][j][5]; // Ajusta direção 6
    }

    // Condição de contorno inferior, adiabático, bounce-back (s distribuições nas direções que apontam para fora do domínio são definidas como as mesmas das células logo acima, simulando uma reflexão perfeita)
    for (int i = 0; i < M; i++) {
        f[i][0][1] = f[i][1][1]; // Reflete as direções para a parte inferior
        f[i][0][4] = f[i][1][4]; // Reflete a direção 5
        f[i][0][5] = f[i][1][5]; // Reflete a direção 6
    }

    // Condição de contorno superior, T=0.0
    for (int i = 0; i < M; i++) {
        f[i][N-1][6] = -f[i][N-1][4]; // Condição para direção 7
        f[i][N-1][1] = -f[i][N-1][3]; // Condição para direção 2
        f[i][N-1][7] = -f[i][N-1][5]; // Condição para direção 8
    }

    // Condição de contorno direita
    for (int j = 0; j < N; j++) {
        f[M-1][j][2] = -f[M-1][j][0]; // Condição para direção 1
        f[M-1][j][6] = -f[M-1][j][4]; // Condição para direção 5
        f[M-1][j][5] = -f[M-1][j][7]; // Condição para direção 6
    }

    // Cálculo da densidade
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            double sumk = 0.0; // Inicializa a soma das direções
            for (int k = 0; k < Q; k++) {
                sumk += f[i][j][k]; // Soma as distribuições
            }
            rho[i][j] = f0[i][j] + sumk; // Calcula a densidade total
        }
    }

    // Rotação da matriz para plotagem de contorno
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            Z[j][i] = rho[i][j]; // Transpõe a matriz de densidade para Z
        }
    }

    // Cálculo de Tm para plotagem
    for (int i = 0; i < N; i++) {
        Tm[i] = rho[i][(N - 1) / 2]; // Seleciona a densidade na linha central
    }

    // Exibir resultados (substitua pela sua própria implementação de gráficos)
    printf("X         T         Flux        Fluxq\n"); // Cabeçalho da tabela
    for (int i = 0; i < M; i++) {
        printf("%lf  %lf  %lf  %lf\n", x[i], Tm[i], flux[i], fluxq[i]); // Exibe as coordenadas x e Tm
    }

    // Salvar os resultados em um arquivo
    FILE *fp;
    fp = fopen("resultado.dat", "w"); // Abre um arquivo para escrita

    if (fp == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        return 1; // Retorna erro se não conseguir abrir
    }

    // Escreve os dados de x e Tm no arquivo
    for (int i = 0; i < M; i++) {
        fprintf(fp, "%lf %lf\n", x[i], Tm[i]); // Formato: x T
    }
    
    fclose(fp); // Fecha o arquivo

    return 0; // Fim da função principal
}
