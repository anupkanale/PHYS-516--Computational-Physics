// Program to perform Monte Carlo Simulations of stock-price using the
// Black-Scholes equation. Based on the random walk example as in diffuse.c

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

double rand_normal();

int main(){
  int N; /* Number of walkers */
  int day,investor,k;
  int hist[60];
	double mu = 0.14; // expected return p.a in %
	double dt=0.00274;
	double sigma=0.2; // Standard deviation of returns
	double stockPrice;

	FILE *fp1, *fp2;
	char filename[50];
	char prefix[17] = "stockData_walker";
	char suffix[4] = ".txt";
	char strN[5];
	
  /* Input parameters */
  printf("Input the number of walkers\n");
  scanf("%d",&N);

  for (k=0; k<60; k++)
    hist[k] = 0;


	for (investor=1; investor<=N; investor++) {
		sprintf(filename, "%s%d%s", prefix,investor,suffix);
		fp1 = fopen(filename, "w");

		stockPrice = 20;
		for (day=1; day<365; day++){
			stockPrice += stockPrice*(mu*dt + sigma* sqrt(dt)* rand_normal());
			fprintf(fp1, "%d %f  \n", day, stockPrice);
		}
		k = round(stockPrice);
		hist[k] += 1;
		fclose(fp1);
	}


	fp2 = fopen("histogramData.txt", "w");
  for (k=0; k<60; k++)
		fprintf(fp2, "%d %d \n", k, hist[k]);
	fclose(fp2);
}

double rand_normal() {
	double r1,r2,eps;
	r1 = rand()/(double) RAND_MAX;
	r2 = rand()/(double) RAND_MAX;
	eps = sqrt(-2* log(r1))*cos(2*M_PI*r2);
	return eps;
}
