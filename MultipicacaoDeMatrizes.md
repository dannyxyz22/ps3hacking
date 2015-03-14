# Objetivo #

Verificar o ganho de tempo do processamento paralelo nas SPEs do Cell em relação ao processamento apenas num PPE. Será realizada a multiplicação de duas matrizes 100000 vezes, para que mesmo usando matrizes pequenas, o ganho de tempo de execução possa ser observado.


# Códigos e Comentários #

Serão utilizadas matrizes 8x8, para que cada SPE seja responsável pelo resultado de uma linha da matriz final.

O código abaixo realiza a multiplicação de matrizes da forma convencional, apenas PPE, e mostra o tempo de execução.

## Programa PPE - serial ##
```
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define SPE_THREADS 8
 
typedef struct
{
    int vezes, n, p, num;
    float linmat1[SPE_THREADS];
    float mat2[SPE_THREADS][SPE_THREADS];
    float linresult[SPE_THREADS];
} __attribute__((aligned(128))) context;
 
int main(
{
    clock_t start, end;
    double elapsed;
 
    context ctxs[SPE_THREADS];
    int i, j, k, l, m, n, p, vezes=100000;
 
    scanf("%d %d %d", &m, &n, &p); //mat1(m x n) mat2(n x p) mat3(m x p)
    m = n = p = SPE_THREADS;
 
    //le mat1
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            scanf("%f", &ctxs[i].linmat1[j]);
 
    //le mat2
    for(j=0; j<n; j++)
        for(k=0; k<p; k++)
            scanf("%f", &ctxs[0].mat2[j][k]);
 
    start = clock();
    //processamento
    for(l=0; l<vezes;l++)
    {
        for(i=0; i<m; i++)
            for(j=0; j<p; j++)
                ctxs[i].linresult[j]=0;
 
        for(i=0; i<m; i++)
            for(j=0; j<n; j++)
                for(k=0; k<p; k++)
                    ctxs[i].linresult[k] += ctxs[i].linmat1[j]*ctxs[0].mat2[j][k];
    }
 
    //imprime a matriz resultado
    for(i=0; i<m; i++, printf("\n"))
        for(j=0; j<p; j++)
            printf("%.2f ", ctxs[i].linresult[j]);
 
    //calcula e imprime o tempo de execução
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Tempo de execucao: %f\n",elapsed);
 
    return 0;
}
```

A estrutura
```
typedef struct
{
    int vezes, n, p, num;
    float linmat1[SPE_THREADS];
    float mat2[SPE_THREADS][SPE_THREADS];
    float linresult[SPE_THREADS];
} __attribute__((aligned(128))) context;
```
será utilizada tanto no processamento PPE quanto SPE, mas a variável _vezes_ e o comando _attribute__((aligned(128))) context_ serão usados apenas no SPE.

A seguir serão mostrados os programas PPE, que distribuirá o “trabalho” entre as SPEs, e SPE, que realiza o cálculo propriamente dito.

## Programa PPE ##
```
1   #include <stdio.h>
2   #include <stdlib.h>
3   #include <time.h>
4   #include <errno.h>
5   #include <libspe.h>
 
6   extern spe_program_handle_t matrix_spu;
 
7   typedef struct
8   {
9       int vezes, n, p, num;
10      float linmat1[SPE_THREADS];
11      float mat2[SPE_THREADS][SPE_THREADS];
12      float linresult[SPE_THREADS];
13  } __attribute__((aligned(128))) context;
 
14  int main()
15  {
16       int i, j, k, status, m, n, p, vezes=1;
17       clock_t start, end;
18       double elapsed;
19       start = clock();  //inicia a contagem do tempo
 
20       scanf("%d %d %d", &m, &n, &p); //mat1(m x n) mat2(n x p) mat3(m x p)
21       m = n = p = SPE_THREADS; //deixamos o tamanho fixo em 8x8
  
22       context ctxs[m]; //cria o vetor de "contextos"
 
23       //le a matriz1, com uma linha para cada "contexto"
24       for(i=0; i<m; i++)
25           for(j=0; j<n; j++)
26       scanf("%f", &ctxs[i].linmat1[j]);
 
27       //le a matriz2, com cada "contexto" possuindo uma copia completa da mesma
28       for(i=0; i<m; i++)
29           for(j=0; j<n; j++)
30               for(k=0; k<p; k++)
31                   if(i==0) scanf("%f", &ctxs[0].mat2[j][k]);
32                   else ctxs[i].mat2[j][k]=ctxs[0].mat2[j][k];
 
33       //inicia as threads para cada SPE, enviando os dados necessarios
34       speid_t spe_ids[m];
35       for(i=0; i<m; i++)
36       {
37           ctxs[i].num = i;
38           ctxs[i].vezes = vezes;
39           ctxs[i].n = n;
40           ctxs[i].p = p;
41           spe_ids[i] = spe_create_thread(0, &matrix_spu, &ctxs[i], NULL, -1, 0);
42           if (spe_ids[i] == 0)
43           {
44               fprintf(stderr, "Failed spe_create_thread(rc=%d, errno=%d)\n",
45               spe_ids[i], errno);
46               exit(1);
47           }
48       }
 
49       // espera pela execucao completa de cada SPU-thread
50       for (i=0; i<m; i++)
51           (void)spe_wait(spe_ids[i], &status, 0);
  
52       //imprime a matriz resultado
53       for(i=0; i<m; i++, printf("\n"))
54           for(j=0; j<p; j++)
55               printf("%f ", ctxs[i].linresult[j]);
  
56       //calcula e imprime o tempo de execucao
57       end = clock();
58       elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
59       printf("Tempo de execucao: %f\n",elapsed);
60       return 0;
61  }
```

Na linha 34, o tipo `speid_t` é um ponteiro para o Local Store de um SPE. Através dele, pode-se criar uma thread usando o comando `spe_create_thread( )`.

Na linha 51, o comando `(void)spe_wait(spe_ids[i], &status, 0)` espera o término de cada uma das _m_ threads

## Programa SPE ##
```
1   #include <stdio.h>
2   #include <spu_intrinsics.h>
3   #include <cbe_mfc.h>
4   #include<libspe.h>
5   #include<libspe2.h>
 
6   #define SPE_THREADS 8
 
7   typedef struct
8   {
9       int vezes, n, p, num;
10       float linmat1[SPE_THREADS];
11      float mat2[SPE_THREADS][SPE_THREADS];
12      float linresult[SPE_THREADS];
13  } __attribute__((aligned(128))) context;
 
14  int main(unsigned long long id, unsigned long long argp)
15  {
16      int i, j, k;
17      context ctx;
 
18      //inicia o recebimento do "contexto" vindo da memória principal para a LS desta SPU
19      spu_mfcdma32(&ctx, argp, sizeof(context), 1, MFC_GET_CMD);
 
20      //espera o final da transferência DMA
21      spu_writech(MFC_WrTagMask, 1<<1); //1<<1: a transferencia 1 está na consulta
22      spu_mfcstat(MFC_TAG_UPDATE_ALL);  //libera o processador quando a transf. 1 terminar
 
23      //realiza o processamento
24      for(i = 0; i < ctx.vezes; i++)
25      {
26          for(k = 0; k < ctx.p; k++) ctx.linresult[k]=0;
27              for(j = 0; j < ctx.n; j++)
28                  for(k = 0; k < ctx.p; k++)
29                      ctx.linresult[k] += ctx.linmat1[j] * ctx.mat2[j][k];
30      }
 
31      //inicia o envio da "linresult" do LS desta SPU de volta para a memória principal
32      spu_mfcdma32(&ctx.linresult, argp + (4 * sizeof(ctx.num)) + sizeof(ctx.linmat1) + sizeof(ctx.mat2), sizeof(ctx.linresult), 2, MFC_PUT_CMD);
33      spu_writech(MFC_WrTagMask, 1<<2);
34      spu_mfcstat(MFC_TAG_UPDATE_ALL); //espera o final da transferência DMA
 
35      return 0;
36      }
```

A função principal _main_ recebe dois parâmetros: _id_ refere-se ao número identificador desta thread; _argp_ representa o ponteiro, pelo qual foram enviados, na chamada de criação da thread no código PPU, os dados a serem processados.

Na linha 19, o recebimento dos dados enviados na criação da thread é iniciado.

As linhas 21 e 22 garantem que a transferência de dados foi concluída.

O loop da linha 24 é responsável por repetir o cálculo 100000 vezes (_ctx.vezes = 100000_), para tornar o processamento mais demorado e ser possível a observação do tempo de execução.

Nas linhas 32, 33 e 34, realiza-se a transferência da linha de resultados _linresult_ do SPE para o PPE, que imprimira a matriz resultado após receber os dados de cada SPE. A impressão não é realizada diretamente por cada SPE porque nada garante que o processamento em cada SPE terminará na mesma ordem em que foi iniciado, ou seja, a ordem das linhas da matriz resultado poderia estar trocada.

# Considerações Finais #

A sintaxe das funções utilizadas e um maior detalhamento do código podem ser encontrados no material usado como bibliografia e base deste post.

Não foi conseguida a compilação do programa SPE, pois algumas das bibliotecas necessárias (_spu\_intrinsics.h_ e _cbe\_mfc.h_) não estão na versão do SDK utilizada ou não estão numa pasta em que o compilador consiga visualizar. Futuramente este problema será solucionado e o resultado do tempo de execução de cada um dos programas, serial e paralelo, será comparado.

# Bibliografia #

FIALHO A. R. S.; MORAIS F.; PLIEUTAUD J. P.; EIZAK L.; Cell Broadband Engine Tutorial. [Online; acessado em 20 de junho de 2008. Disponível em [link](http://www.pad.lsi.usp.br/joomla/index.php?option=com_docman&task=doc_download&gid=79&Itemid=89)]

IBM. Cell Broadband Engine Resource Center [Online; acessado em 20 de junho de 2008. Disponível em [link](http://www-128.ibm.com/developerworks/power/cell/documents.html?S_TACT=105AGX16&S_CMP=LP)]