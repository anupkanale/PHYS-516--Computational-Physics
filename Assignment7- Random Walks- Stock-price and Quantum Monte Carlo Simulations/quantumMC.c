#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// define function prototypes
void intialize();
void walk();
void data();

// Constants
int N0=50;	 		//initial walkers
double ds=0.1;  //random walk step
double dt=0.01; //time step (=ds^2)
int mcs=500; 		//MC steps

// Variables
int N;			 //# of walkers
double Vref; //ref energy

int main() {

	initialize(N);
	for(step=1; step<=0.5*mcs; step++)
		walk();		 // thermalization

	for (step=1; step<=mcs; step++) {
		walk();
		data();
	}
	//plot_psi();

  return 0;
}

// Function definions

void initilize(){
	N=N0
	for (i=1; i<N; i++)
		x[i] = uniform random number # in [-1,1]

	V_ref ← 1/N \sigma_i=1^N V(x_i) // = <V>, this is the harmonic potential x_i^2/2
}

void walk(){
Nin ← N;
Vsum ← 0;
for (i=Nin; i>=1; i--){ //Going ulta!!
	// Random Walk as in diffuse.c
	if(rand()%2==0) x[i] += ds;
	else x[i] -= ds;
	// Birth or death
	dV = V(x_i) – Vref; //here V(x_i) = x_i^2/2
	
	if (dV<0.0)
		if ( rand()/(double)RANDMAX < -dV*dt)
			++N;
			x[N] = x[i];
			Vsum += 2*V(x_i);
		else
			Vsum += V(x_i);
	else
		if ( rand()/(double)RANDMAX < dV*dt)
			x[i] = x[N];
			N--;
		else
			Vsum += V(x_i);
end for

Vavg = Vsum/N;
Vref = Vavg – (N-N0)/(N0*dt);
}

void data(){
	for(i=1; i<=N; i++)
		++psi[(int) (x_i + xshift + 0.5*binsize)/binsize];
}
