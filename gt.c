#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <locale.h>
#define T1 0.000009
#define T2 0.0001
#define T3 0.0002
#define ERR 1

///--Estruturas de armazenamento ------------------------------------------
typedef struct {
	double pot_min, pot_max;
	double a, b, c;
	double fator_p, err;

}u_geradora;

///--Cabeçalho de Funções -------------------------------------------------
bool calc_erro (u_geradora geradoras[], int qt_u);
bool calc_min (u_geradora geradoras[], int qt_u);
double perda_transmissao(double p1, double p2, double p3);
double newton_raphson (u_geradora u, double pn);
double fator_participacao (u_geradora u, double so_s);
double serie_soma (u_geradora geradoras[], int qt_u);
double f_n (u_geradora u, double p);
double df_n (u_geradora u, double p);
double ddf_n (u_geradora u);

int main(){
	
	//Setando a linguagem em Português
	setlocale (LC_ALL, "Portuguese");
	
	double carga_t;
	int qt_u;
	
	printf ("Insira a quantidade de Geradoras Térmicas: \n");
	scanf ("%d", &qt_u);
	
	printf ("Insira carga em MW a ser suprida: \n");
	scanf ("%lf", &carga_t);
	
	//Criação do vetor de unidades geradoras
	u_geradora geradoras[qt_u];
	short x;
	
	//Captura dos valores
	for (x = 0; x < qt_u; x++) {
		fflush (stdin);
		printf ("Insira a potência mínima da Unidade Geradora %d: \n", x+1);
		scanf ("%lf", &geradoras[x].pot_min);
		
		printf ("Insira a potência máxima da Unidade Geradora %d: \n", x+1);
		scanf ("%lf", &geradoras[x].pot_max);
		
		printf ("Insira 'a', 'b' e 'c' da equação de Fator de potência da Unidade Geradora %d, respectivamente: \n", x+1);
		
		scanf ("%lf", &geradoras[x].a); fflush (stdin);
		scanf ("%lf", &geradoras[x].b); fflush (stdin);
		scanf ("%lf", &geradoras[x].c); fflush (stdin);
		
	}
	
	//Cálculo da Série 
	double so_s = serie_soma (geradoras, qt_u);
	
	//Calculo do Fator de participação para cada unidade geradora
	for (x = 0; x < qt_u; x++) {
		double fp = fator_participacao (geradoras[x], so_s);
		printf("%lf  -- ", fp);
		double r = carga_t * fp;
		geradoras[x].fator_p = r;
		geradoras[x].err = 10;
		printf("%lf\n", geradoras[x].fator_p);
	}
	printf("\n\n");
	
	short it = 1;
	//Implementação do método de Newton-Raphson
	while (calc_erro (geradoras, qt_u) && calc_min (geradoras, qt_u)) {
		
		printf("\n\nIteração %d\n\n", it);
		short z;
		for (z = 0; z < qt_u; z++) {
			double p_novo = newton_raphson (geradoras[z], geradoras[z].fator_p);
			geradoras[z].err = fabs (p_novo - geradoras[z].fator_p);
			geradoras[z].fator_p = p_novo;
			printf("Geradora %d == %lf\n", z, geradoras[z].fator_p);
			printf("Erro %d == %lf\n", z, geradoras[z].err);			
		}
		it++;
	}
	
	if (qt_u == 3) {	
		double pt = perda_transmissao (geradoras[0].fator_p, geradoras[1].fator_p, geradoras[2].fator_p);
		printf("\n\nPerdas Totais de: %lf MW", pt);
	}
	return 0;
}

///--Implementação das Funções -------------------------------------------
bool calc_erro (u_geradora geradoras[], int qt_u) {
	short x;
	for (x = 0; x < qt_u; x++) 
		if (geradoras[x].err > ERR) return true;
	return false;
}

bool calc_min (u_geradora geradoras[], int qt_u) {
	short x;
	for (x = 0; x < qt_u; x++)
		if (geradoras[x].fator_p <= geradoras[x].pot_min) return false;
	return true;
}

double perda_transmissao(double p1, double p2, double p3) {
	return T1 * pow(p1, 2) + T2 * pow(p2, 2) + T3 * pow(p3, 2);
}

double newton_raphson (u_geradora u, double pn) {
	return fabs (pn - ( f_n (u, pn) / df_n (u, pn) ));
}

double fator_participacao (u_geradora u, double so_s) {
	return ((1 / ddf_n (u)) / so_s ); 
}

double serie_soma (u_geradora geradoras[], int qt_u) {
	
	double res = 0;
	short x;
	for (x = 0; x < qt_u; x++) 
		res += (1 / ddf_n (geradoras[x]));	
	
	return res;
}

double f_n (u_geradora u, double p) {
	return pow(p, 2) * u.a + p * u.b + u.c;
}

double df_n (u_geradora u, double p) {
	return p * 2 * u.a + u.b;
}

double ddf_n (u_geradora u) {
	return 2 * u.a;
}
