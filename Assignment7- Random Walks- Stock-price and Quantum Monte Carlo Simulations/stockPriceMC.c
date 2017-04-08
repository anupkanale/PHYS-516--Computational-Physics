// Program to perform Monte Carlo Simulations of stock-price using the
// Black-Scholes equation

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double rand_normal();

int main(){
  int N; /* Number of walkers */
  int day,investor,k;
  int hist[1001];

	double mu = 0.14, dt=0.00274, sigma=0.2;
	double stockPrice;

  /* Input parameters */
  printf("Input the number of walkers\n");
  scanf("%d",&N);

  for (k=0; k<=1000; k++)
    hist[k] = 0.0;

	FILE *fp;
	fp = fopen("stockData.txt", "w");
	for (investor=1; investor<=N; investor++) {
	stockPrice = 20;
		for (day=1; day<365; day++){
			stockPrice += stockPrice*(mu*dt + sigma* sqrt(dt)* rand_normal());
			fprintf(fp, "%d %f  \n", day, stockPrice);
		}
	}

	fclose(fp);
}

double rand_normal() {
	double r1,r2,eps;
	r1 = rand()/(double) RAND_MAX;
	r2 = rand()/(double) RAND_MAX;
	eps = sqrt(-2* log(r1))*cos(2*M_PI*r2);
	return eps;
}
