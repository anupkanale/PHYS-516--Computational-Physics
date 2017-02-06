/* Monte Carlo Simulation of 2D Ising Model */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Define global variables
#define L 20 //lattice size
int s[L][L]; //Spins s[i][j]=+-1
double exp_dV[2][5];
double JdivT; // J/kBT
double HdivT; // H/kBT

void table_set(){ //function to compute the exponential term
  int sDash,sneighbor,k,l;
  for (k=0;k<2;k++){
  	sDash = 2*k-1;
    for(l=0;l<5;l++){
    	sneighbor = 2*l-4;
    	exp_dV[k][l] = exp(2*sDash*(JdivT*sneighbor + HdivT)); // check formula 
    }
  }
}

int main() {
	double runM;
	double sumM=0.0, sumM2=0.0;
	double exp_val, avgM, sigM;
	int snew,sneighbor,Sta_step;
	int i,j,step,k,l,im,ip,jm,jp;

	printf("Input J/kBT\n");
	scanf("%le",&JdivT);

	HdivT = 0.0;
	Sta_step = 2000000;

	table_set();	// Set up the look-up table for the exponent calculation

	for(i=0;i<5;i++) {			//Cold start- start with all spins up configuration
		for(j=0;j<5;j++) {
			s[i][j] = 1;
		}
	}
	runM=1.0*L*L;

	for(step=0; step<Sta_step; step++) {
		i=rand()%L;
		j=rand()%L;
		snew=-s[i][j];

		// Figure out which element of the table is to be looked up
		im = (i + L - 1) % L;
		ip = (i + 1) % L;
		jm = (j + L - 1) % L;
		jp = (j + 1) % L;

		k = (snew+1)/2;
		sneighbor = s[im][j] + s[ip][j] + s[i][jm] + s[i][jp];
		l = (sneighbor+4)/2;

		//Change in Pot Energy wth flip
		exp_val=exp_dV[k][l];

		// Accept of reject flip
		if (exp_val>1.0) {
			s[i][j] = snew;
			runM += 2*snew; //update value of magnetization
		}
		else if(rand()/(double)RAND_MAX < exp_val) {
			s[i][j] = snew;
			runM += 2*snew; //update value of magnetization
		}

		sumM += runM;
		sumM2 += runM*runM;
	}
	avgM = sumM/Sta_step;					//Mean Magnetization
	sigM = sqrt(sumM2/Sta_step-avgM*avgM);	//Standard deviation

	printf("Mean Magnetization %le \n", avgM);
	printf("Standard Deviation %le \n", sigM);
}