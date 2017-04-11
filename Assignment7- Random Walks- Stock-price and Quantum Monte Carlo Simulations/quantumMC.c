/*
Quantum Monte Carlo Algorithm
--------------------------------------
1. Place N0 walkers at the initial set of positions xi.
2. Compute the reference energy, Vref = Σi V(xi)/N0.
3. For each walker,
	a) Randomly move the walker to the right or left by a fixed step length ∆s.
	b) Compute ΔV = V(x) − Vref and a random number r in the range [0, 1]. If ΔV > 0 and r <
	ΔVΔτ, then remove the walker. If ΔV <0 and r < −ΔVΔτ, then add another walker at x.
	Otherwise, just leave the walker at x.
4. Compute the mean potential energy (6) and the actual number of random walkers. The new
reference potential is given by Eq. (7). The average 〈V〉 is an estimate of the ground state
energy.
5. Repeat steps 3−4 until the estimates of the ground state energy 〈V〉 have reached a steady
state value with only random fluctuations. Average 〈V〉 over many Monte Carlo steps to
compute the ground state energy.
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// define function prototypes
void initialize();
void walk();
void data();

// Constants
#define N0 50	 	//initial walkers
#define mcs 500	//MC steps
#define MAXX 100
#define maxpsi 100

double ds=0.1;  //random walk step
double dt=0.01; //time step (=ds^2)


// Variables
int N;		 //# of walkers
double Vref; //ref energy

double x[MAXX+1];// position of each random walker
double psi[maxpsi+1]; // histogram

int main() {
	int step;
	FILE *fp;

  srand((unsigned)time((long *)0)); /* Initialize the random-number sequence */

	initialize();
	printf("%f \n", Vref);
	for(step=1; step<=0.4*mcs; step++){
		walk();		 // thermalization
	}

	fp = fopen("vRef.txt", "w");
	for (step=1; step<=mcs; step++) {
		walk();
		fprintf(fp, "%d %f \n",step, Vref);
		data();
	}
	fclose(fp);
	//plot_psi();

 	return 0;
}

// Function definitions
void initialize() {
	int i;
	N = N0;

	// Assign position to N0 random walkers
	for (i=0; i<N; i++)
		x[i] = -1+2*rand()/(double) RAND_MAX ; //uniform random number # in [-1,1]

	// Calculate reference energy
	for (i=0; i<N; i++)
		Vref += pow(x[i],2)/2; // = <V>, this is the harmonic potential x_i^2/2
	Vref = Vref/N;
}

void walk() {
	int i;
	int Nin = N;
	double Vavg, dV, Vsum=0.0;

	for (i=Nin; i>=1; i--) { //Going ulta!!
	// Random Walk as in diffuse.c
		if(rand()%2==0) x[i] += ds;
		else x[i] -= ds;

		// Birth or death
		dV = pow(x[i],2)/2 - Vref;
		if (dV<0.0)
			if (rand()/(double)RAND_MAX < -dV*dt){
				N++;
				x[N] = x[i];
				Vsum += 2*pow(x[i],2)/2;} //here V(x_i) = x_i^2/2}
			else
				Vsum += pow(x[i],2)/2;
		else
			if ( rand()/(double)RAND_MAX < dV*dt){
				x[i] = x[N];
				N--;}
			else
				Vsum += pow(x[i],2)/2;
	}

	Vavg = Vsum/N;
	Vref = Vavg - (N-N0)/(double)(N0*dt);
}

void data(){
	double xshift, binsize;
	xshift = binsize*maxpsi*0.5;
	binsize = 2*ds;
	for(int i=1; i<=N; i++){}
		//++psi[(int) (x[i] + xshift + 0.5*binsize)/binsize];
}
