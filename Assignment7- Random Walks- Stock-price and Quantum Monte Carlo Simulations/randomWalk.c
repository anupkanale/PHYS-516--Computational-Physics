// Program to perform Monte Carlo Simulations of stock-price using the
// Black-Scholes equation

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main(){
	int maxStep = 500, N_walker = 1;
	int position;


	FILE *fp;
	fp = fopen("drunkGuy.txt", "w");
	fprintf(fp, "Positions of the drunk guy\n");
	fprintf(fp, "/*--------------------------------------*/\n");
	for (int walker=0; walker<N_walker; walker++) {
		position=0;
	    for (int step=0; step<maxStep; step++) {
	        if ((double) rand()/RAND_MAX>0.5)	position += 1;
	        else						position -= 1;

	        fprintf(fp, "%d  \n", position);
   	    }
	}
	fclose(fp);
}