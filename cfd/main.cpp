#include <stdio.h>
#include <vector>
#include <malloc.h>
#include <cmath>
#include <rapidmind/platform.hpp>
#include <rapidmind/shortcuts.hpp>
using namespace rapidmind;

using namespace std;

class No{
public:
	double x,y,phi,delta,velocidade;

	No(double x, double y){
		this->x = x;
		this->y = y;
	}
};


enum Metodo { Jacobi, GaussSeidel, SOR, LineGaussSeidel, SLOR };

Timer timerStart, timerEnd;

//estrutura utilizada para armazenar uma malha 2d de nos
//o construtor espera uma equacao da forma a x^2 + b x + c = 0
//para que todos os pontos menores que estes nao sejam incluidos na malha

class Malha{
private:
	int ILE,ITE,IMAX,JMAX;
	double XSF,YSF, t, Uinf, eps, relaxacaoSor;
    int count;

	vector< vector <No> > nos;
	vector< No > contorno;
	vector< double > historicoResiduo;//em log na base 10
	vector< double > cpAerofolio; //calculo do cP apenas no aerofolio

	void construirMalha(){

		nos = vector< vector <No> >(IMAX+1, vector<No>(JMAX+1, No(0.0,0.0) )) ;
		contorno = vector< No >(IMAX+1,No(0.0,0.0));
		cpAerofolio = vector < double > (ITE - ILE +1, 0.0);

        for(int i=1; i<=IMAX;i++){
			for(int j=1;j<=JMAX;j++){
                printf("%lf  %lf\n", nos[i][j].x, nos[i][j].y);
            }
		}

		double deltaX = 1.0/ (ITE-ILE);

		printf("%lf %d %d\n", deltaX, ITE, ILE);

		printf("%d  ", ITE);
        printf("%d\n", ILE);

		//note que deve-se primeiro preencher de ILE a ITE
		//pois os outros preenchimentos sao dependentes deste
		for(int i=ILE;i<=ITE;i++){
			for(int j=1;j<=JMAX;j++){
				nos[i][j].x = (i-ILE)* deltaX;
				printf("%lf\n", nos[i][j].x);
			}
		}
		printf("%lf %d %d", deltaX, ITE, ILE);
		for(int i=ITE+1;i<=IMAX;i++){
			for(int j=1;j<=JMAX;j++){
				nos[i][j].x = nos[i-1][j].x + ( nos[i-1][j].x - nos[i-2][j].x)*XSF;
				printf("[%d][%d] %lf\n", i, j, nos[i][j].x);
			}
		}

		//tomar o cuidado de preencher da direita para a esquerda
		for(int i=ILE-1;i>=1;i--){
			for(int j=1;j<=JMAX;j++){
				nos[i][j].x = nos[i+1][j].x + ( nos[i+1][j].x - nos[i+2][j].x)*XSF;
			}
		}

		//malha em y
		for(int i=1; i<=IMAX; i++){
			nos[i][1].y= -deltaX/2.0;
			nos[i][2].y=  deltaX/2.0;
			printf("%lf %lf\n", nos[i][1].y, nos[i][2].y);
		}
		for(int i=1; i<=IMAX;i++){
			for(int j=3;j<=JMAX;j++){
				nos[i][j].y= nos[i][j-1].y + (nos[i][j-1].y- nos[i][j-2].y)*YSF;
				printf("[%d][%d]  %lf\n", i, j, nos[i][j].y);
			}
		}

		imprimirMalha();
		imprimirAerofolio();

	}
	//as condicoes iniciais serao da forma phiInf = Uinf*x

	void aplicarCondicoesIniciais(){
	    printf("WPHI\n");
		for(int i=1;i<=IMAX; i++){
			for(int j=1;j<=JMAX;j++){
				nos[i][j].phi = Uinf*nos[i][j].x;
				printf("[%d][%d]  %lf\n", i, j, nos[i][j].phi);
			}
		}
        printf("CONTORNO\n");
		for(int i=1;i<=IMAX;i++){
			//para j=1.5, temos informacao sobre dphidy
			if((i<ILE)||(i>ITE)){
				//condicao de simetria do escoamento para alpha=0
				contorno[i].phi=0.0;
				printf("%d  %lf\n", i, contorno[i].phi);
			}
			else{
				//no perfil/corda temos que
				//dphidy = Uinf * dy/dx
				//dy/dx = 2t - 4tx
				//note que poderiamos pegar qualquer j para nos[i][1].x
				//visto que x eh constante na direcao y
				contorno[i].phi=Uinf* (2*t - 4*t*nos[i][1].x);
                printf("%d  %lf\n", i, contorno[i].phi);
			}
		}
	}

	void aplicarDeltas(){
		for(int j=2;j<JMAX;j++){
			for(int i=2;i<IMAX;i++){
				nos[i][j].phi += nos[i][j].delta;
			}
		}
	}
	//os metodos Jacobi, GaussSeidel e SOR serao iterados na direcao y
	bool iterarDirecaoY(Metodo metodo){
		bool terminou = true;

		//observar nao iterar no contorno
		double resMax = 0;

		for(int j=2;j<JMAX;j++){
			for(int i=2;i<IMAX;i++){
				double N;
				//conforme definido na handout4
				double deltaX = (nos[i+1][j].x - nos[i-1][j].x)/2.0;
				double deltaY = (nos[i][j+1].y - nos[i][j-1].y)/2.0;

				double Lij = 2.0/(nos[i+1][j].x-nos[i-1][j].x)*
				((nos[i+1][j].phi-nos [i] [j].phi)/(nos[i+1][j].x-nos[i][j].x)-
						((nos[i]  [j].phi-nos[i-1][j].phi)/(nos[i][j].x-nos[i-1][j].x))
				)+
				2.0/(nos[i][j+1].y-nos[i][j-1].y)*
				( (nos[i][j+1].phi-nos[i][j].phi)/(nos[i][j+1].y-nos[i][j].y)-
						((nos[i][j].phi - nos[i][j-1].phi)/(nos[i][j].y-nos[i][j-1].y))
				);

				if(fabs(Lij)>resMax){
					resMax = fabs(Lij);
				}
				//printf("%+.9lf ",Lij);
				switch(metodo){
				case Jacobi:
					N = -2.0/(deltaX*deltaX) -2.0/(deltaY*deltaY);
					nos[i][j].delta = -Lij/N;
					break;
				case GaussSeidel:
					N = -2.0/(deltaX*deltaX) -2.0/(deltaY*deltaY);
					nos[i][j].delta = (-Lij-nos[i-1][j].delta-nos[i][j-1].delta)/N;
					break;
				case SOR:
					double r = relaxacaoSor;
					N = -(2.0/r)/(deltaX*deltaX) -(2.0/r)/(deltaY*deltaY);
					nos[i][j].delta = (-Lij-nos[i-1][j].delta-nos[i][j-1].delta)/N;
					break;
				}

				if((terminou)&&(fabs(Lij)>=1.0e-8)){
					//	printf("Ainda nao (%d,%d) %.30lf\n",i,j,fabs(Lij/N));
					terminou = false;
				}

			}
		}

		aplicarDeltas();

		double log10Res = log(resMax)/log(10);
		historicoResiduo.push_back(log10Res);
		//printf("\n");
		return terminou;
	}

	bool iterarDirecaoX(Metodo metodo){
		// malha(a,b,c,d) indiceMalha
		// contorno(n+2)   JMAX
		// n+1             ...
		// n               ...
		// ...             ...
		// ...             ...
		// 2                2
		// contorno(1)      1

		double resMax = 0;

		bool terminou = true;

		vector <double>  a = vector< double >(JMAX+1, 0.0) ;
		vector <double>  b = vector< double >(JMAX+1, 0.0) ;
		vector <double>  c = vector< double >(JMAX+1, 0.0) ;
		vector <double>  d = vector< double >(JMAX+1, 0.0) ;
		vector <double>  as = vector< double >(JMAX-1, 0.0) ;
		vector <double>  bs = vector< double >(JMAX-1, 0.0) ;
		vector <double>  cs = vector< double >(JMAX-1, 0.0) ;
		vector <double>  ds = vector< double >(JMAX-1, 0.0) ;
		vector <double>  x = vector< double >(JMAX-1, 0.0) ;


		for(int i=2;i<IMAX;i++){
			//for(int j=2;j<JMAX;j++){

			int n  = JMAX-2;//iteracao vai para j de 2 a n+1(JMAX-1), para i de 2 a (IMAX-1)

			double Li2,Lin1,deltaX,deltaY;
			for(int j = 2;j<= JMAX-1;j++){
				double Lij = 2.0/(nos[i+1][j].x-nos[i-1][j].x)*
				((nos[i+1][j].phi-nos [i] [j].phi)/(nos[i+1][j].x-nos[i][j].x)-
						((nos[i]  [j].phi-nos[i-1][j].phi)/(nos[i][j].x-nos[i-1][j].x))
				)+
				2.0/(nos[i][j+1].y-nos[i][j-1].y)*
				( (nos[i][j+1].phi-nos[i][j].phi)/(nos[i][j+1].y-nos[i][j].y)-
						((nos[i][j].phi - nos[i][j-1].phi)/(nos[i][j].y-nos[i][j-1].y))
				);

				if(fabs(Lij)>resMax){
					resMax = fabs(Lij);
				}

				//conforme definido na handout4
				deltaX = (nos[i+1][j].x - nos[i-1][j].x)/2.0;
				deltaY = (nos[i][j+1].y - nos[i][j-1].y)/2.0;
				double r=0.0;
				if(metodo==LineGaussSeidel)
					r = 1.0;
				else
					r = relaxacaoSor;

				a[j] =  (2/r)/((nos[i][j].y-nos[i][j-1].y)* (nos[i][j+1].y - nos[i][j-1].y));
				b[j] = -(2/r)/(deltaX*deltaX) - (2/r)/((nos[i][j+1].y - nos[i][j].y)*(nos[i][j+1].y - nos[i][j-1].y)) - (2/r)/((nos[i][j].y-nos[i][j-1].y)* (nos[i][j+1].y - nos[i][j-1].y));
				c[j] =  (2/r)/((nos[i][j+1].y - nos[i][j].y)*(nos[i][j+1].y - nos[i][j-1].y));
				d[j] = -Lij - nos[i-1][j].delta/(deltaX*deltaX);
				if(j==2) 	Li2=Lij;
				if(j==n+1) 	Lin1 = Lij;
			}
			//d[2] e d[n+1] sao casos especiais

			d[2]   = -Li2  - nos[i-1][2].delta/(deltaX*deltaX) - a[2]   * nos[i][1].delta;
			d[n+1] = -Lin1 - nos[i-1][n+1].delta/(deltaX*deltaX) - c[n+1] * nos[i][n+2].delta;

			//poderiamos tambem dizer que
			a[2] = 0, c[n+1] = 0;
			for(int w=2;w<=n+1;w++){
				as[w-1]= a[w];
				bs[w-1]= b[w];
				cs[w-1]= c[w];
				ds[w-1]= d[w];
			}

			tridiag(as,bs,cs,ds,x);
			/*
			if(i==2){
				for(int w=1;w<=n;w++)
				printf("delta(j %d) %lf\n",w,x[w]);
			}*/

			for(int w=1;w<=n;w++){
				nos[i][w+1].delta = x[w];
			}
		}

		double log10Res = log(resMax)/log(10);
		historicoResiduo.push_back(log10Res);

		aplicarDeltas();
		return false;//terminou;
	}

	//faz a iteracao na malha conforme metodo
	bool iterarMalha(Metodo metodo){
		bool terminou = true;
		aplicarCondicoesDeContorno();

		switch (metodo){
		case Jacobi:
		case GaussSeidel:
		case SOR:
			terminou = iterarDirecaoY(metodo);
			break;
		case LineGaussSeidel:
		case SLOR:
			terminou = iterarDirecaoX(metodo);
			break;
		}
		return terminou;
	}

	//aplica as condicoes de contorno nas
	//fronteiras de entrada, saida, superior e inferior
	void aplicarCondicoesDeContorno(){
		//a condicao de contorno nas fronteiras de entrada,
		//saida e superior sao iguais as condicoes iniciais,
		//portanto, somente a fronteira inferior sera atualizada

		//contornos esquerdo e direito
		/**ITERAR CONTORNO*/
		printf("ITERAR CONTORNO\n");
		for(int j=1;j<=JMAX;j++){
			nos[1][j].delta    = 0;
			nos[IMAX][j].delta = 0;
		}
		//contorno superior
		for(int i=1;i<=IMAX;i++){
			nos[i][JMAX].delta = 0;
		}

		//Conforme a lista, o perfil esta sendo representado
		//por sua corda e temos que:
		//phi(i,1) = phi(i,2)-(y2-y1)dphidy(i,3/2)
		//lembrando que o y em 3/2 vale zero, pois
		//y(3/2) = (y(1)+y(2))/2 = 0/2 = 0
		//y2-y1 = deltaX
		double deltaX = 1.0/(ITE-ILE);
		printf("DX %lf\n", deltaX);

		/**ITERAR RESTO*/
		printf("ITERAR RESTO %d\n", count);
		for(int i=1;i<=IMAX;i++){
			//lembrando que contorno[i].phi = dphidy, por construcao
			//TODO: confirmar se nos[i][2] eh realmente da iteracao anterior,
			//como esta sendo feito aqui
			double delta = nos[i][2].phi - ((deltaX) * contorno[i].phi) - nos[i][1].phi;
			printf("PHI = %lf   DX = %lf  CONT = %lf  PHI = %lf\n", nos[i][2].phi, deltaX, contorno[i].phi, nos[i][1].phi);
			nos[i][1].delta = delta;
			nos[i][1].phi +=delta;
			printf("%d  %lf  %lf\n", i, nos[i][1].delta, nos[i][1].phi);
			//printf("Mydelta                %lf\n",delta);
			//nos[i][1].phi += delta;
		}
	}
	double getEpsilonMaquina(){
		double eps = 1.0;
		do{
			printf("%.30lf\n",eps);
			eps = eps/2.0;
		}while( (double)(1.0 + eps/2.0)!=1.0);

		return eps;
	}

public:
	Malha(int ILE, int ITE, int IMAX, int JMAX, double XSF, double YSF,
			double t, double Uinf){
		this->ILE = ILE; 	this->ITE = ITE; 	this->IMAX = IMAX;
		this->JMAX = JMAX;	this->XSF = XSF;	this->YSF = YSF;
		this->t = t,		this->Uinf=Uinf;

		construirMalha();


	}

	//assumindo a[1]=0
	void tridiag(vector<double> a, vector<double> b,vector<double> c,vector<double> d,
			vector<double>& x){

		int n = a.size()-1;
		//vector<double> bl = vector<double> (n+1,0.0);
		//vector<double> dl = vector<double> (n+1,0.0);

		for(int k=2;k<=n;k++){
			double m   = a[k]/b[k-1];
			b[k] = b[k] - m*c[k-1];
			d[k] = d[k] - m*d[k-1];
			//printf("bl(%d) %lf dl(%d) %lf\n",k,b[k],k,d[k]);
		}
		x[n] = d[n]/b[n]; //verificar se eh b ou bl
		for(int k=n-1;k>=1;k--){
			x[k]= (d[k] - c[k]*x[k+1])/b[k]; //verificar se eh bk ou bkl
		}
		return;
	}

	void imprimirAerofolio(){
		FILE* faero = fopen("aero.dat","w");
		for(int i=0;i<50;i++){
			double x = i/50.0;
			double a = -2*t, b = 2*t,c=0;
			double y = a*x*x+b*x+c;
			fprintf(faero,"%lf %lf\n",x,y);
		}

		fclose(faero);
	}

	void imprimirMalha(){
		FILE* fmalha = fopen("malha.dat","w");
		for(int j=1;j<=JMAX;j++){
			for(int i=1;i<=IMAX;i++){
				fprintf(fmalha,"%+5.3lf %+5.3lf\n",nos[i][j].x, nos[i][j].y);
			}
		}
		fclose(fmalha);
	}
	void imprimirPotenciais(const char* nomeDoArquivo){
		FILE* fpotencial = fopen(nomeDoArquivo,"w");
		for(int j=1;j<=JMAX;j++){
			for(int i=1;i<=IMAX;i++){
				fprintf(fpotencial,"%+5.5lf %+5.5lf %+5.5lf\n",nos[i][j].x, nos[i][j].y,nos[i][j].phi);
			}
			fprintf(fpotencial,"\n");
		}
		fclose(fpotencial);

	}

	void imprimirResiduo(const char* nomeDoArquivo){
		FILE* fresiduo = fopen(nomeDoArquivo,"w");
		for(int i=0;i<historicoResiduo.size();i++){
			fprintf(fresiduo,"%d %lf\n",i+1,historicoResiduo[i]);
		}
		fclose(fresiduo);

	}

	void imprimirCoeficientesDeVelocidade(const char* nomeDoArquivo){
		FILE* fCoeficientes = fopen(nomeDoArquivo,"w");
		for(int j=2;j<JMAX;j++){
			for(int i=2;i<IMAX;i++){
				fprintf(fCoeficientes,"%+5.5lf %+5.5lf %+5.5lf\n",nos[i][j].x, nos[i][j].y,nos[i][j].velocidade);
			}
			fprintf(fCoeficientes,"\n");
		}

		fclose(fCoeficientes);
	}
	void imprimirCoeficientesDePressao(const char* nomeDoArquivo){
		FILE* fCoeficientes = fopen(nomeDoArquivo,"w");
		for(int j=2;j<JMAX;j++){
			for(int i=2;i<IMAX;i++){
				double cp = 1 - (nos[i][j].velocidade*nos[i][j].velocidade/(Uinf*Uinf));
				fprintf(fCoeficientes,"%+5.5lf %+5.5lf %+5.5lf\n",nos[i][j].x, nos[i][j].y,cp);
			}
			fprintf(fCoeficientes,"\n");
		}

		fclose(fCoeficientes);
	}
	void setRelaxacaoSor(double val){
		relaxacaoSor = val;
	}

	void calcularVelocidades(){
		for(int j=2;j<JMAX;j++){
			for(int i=2;i<IMAX;i++){
				double vx = (nos[i+1][j].phi - nos[i-1][j].phi)/(nos[i+1][j].x - nos[i-1][j].x);
				double vy = (nos[i][j+1].phi - nos[i][j-1].phi)/(nos[i][j+1].y - nos[i][j-1].y);
				nos[i][j].velocidade = sqrt(vx*vx+vy*vy);
			}
		}
	}

	void calcularCPAerofolio (){
		for(int i=ILE;i<=ITE;i++){
			double vxi1 = (nos[i+1][1].phi - nos[i-1][1].phi)/(nos[i+1][1].x - nos[i-1][1].x);
			double vxi2 = (nos[i+1][2].phi - nos[i-1][2].phi)/(nos[i+1][2].x - nos[i-1][2].x);
			double vx = 0.5*(vxi1+vxi2);
			double vy = contorno[i].phi;
			double velocidade = sqrt(vx*vx+vy*vy);
			cpAerofolio[i-ILE]= 1 - (velocidade*velocidade/(Uinf*Uinf));
		}
	}

	void imprimirCpAerofolio(const char* nomeDoArquivo){
		FILE* fcpaer = fopen(nomeDoArquivo,"w");
		for(int i=0;i<cpAerofolio.size();i++){
			fprintf(fcpaer,"%+5.5lf %+5.5lf \n",nos[i+ILE][1].x,-cpAerofolio[i]);
		}

		fclose(fcpaer);
	}

	char* prefixarMetodo(const char* palavra,Metodo metodo){

		char* temp = (char*) malloc(200*sizeof(char));
		switch(metodo){
		case Jacobi:
			sprintf(temp,"Jacobi%s",palavra);
			break;
		case GaussSeidel:
			sprintf(temp,"GaussSeidel%s",palavra);
			break;
		case SOR:
			sprintf(temp,"SOR%s",palavra);
			break;
		case LineGaussSeidel:
			sprintf(temp,"LineGaussSeidel%s",palavra);
			break;
		case SLOR:
			sprintf(temp,"SLOR%s",palavra);
			break;
		}
		return temp;
	}
	void solucionar(Metodo metodo){

		aplicarCondicoesIniciais();

		imprimirPotenciais("potInicial.dat");

		count=0;
		//while(!iterarMalha(metodo) && count < 10000) count++;
		while(count < 10000){
			iterarMalha(metodo);
			count++;
		}

		printf("Iteracoes %d\n",count);

		calcularVelocidades();

		calcularCPAerofolio();

		imprimirCpAerofolio(prefixarMetodo("cpAerofolio.dat",metodo));

		imprimirResiduo(prefixarMetodo("residuo.dat",metodo));

		imprimirPotenciais("potFinal.dat");

		imprimirCoeficientesDeVelocidade("velocidade.dat");

		imprimirCoeficientesDePressao(prefixarMetodo("pressao.dat",metodo));
	}
};

int main(){

	double t = 0.1, Uinf = 1.0;
	double rsor = 0.1;

	finish();
    timerStart = Timer::now();

	Malha m1 = Malha(11,31,41,12,1.25,1.25, t,Uinf);
	m1.setRelaxacaoSor(rsor);
	m1.solucionar(Jacobi);

	finish();
    timerEnd = Timer::now();
    cout << "o tempo de execucao 1 foi de " << (timerEnd - timerStart).milliseconds() << " milissegundos.\n";


	/*Malha m2 = Malha(11,31,41,12,1.25,1.25, t,Uinf);
	m2.setRelaxacaoSor(rsor);
	m2.solucionar(SOR);

	Malha m3 = Malha(11,31,41,12,1.25,1.25, t,Uinf);
	m3.setRelaxacaoSor(rsor);
	m3.solucionar(GaussSeidel);

	Malha m4 = Malha(11,31,41,12,1.25,1.25, t,Uinf);
	m4.setRelaxacaoSor(rsor);
	m4.solucionar(LineGaussSeidel);

	Malha m5 = Malha(11,31,41,12,1.25,1.25, t,Uinf);
	m5.setRelaxacaoSor(rsor);
	m5.solucionar(SLOR);*/

}
