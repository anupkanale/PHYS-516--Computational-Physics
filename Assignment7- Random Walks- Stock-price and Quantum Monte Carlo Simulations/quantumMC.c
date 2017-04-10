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
void intialize();
void walk();
void data();

// Constants
int N0=50;	 	//initial walkers
double ds=0.1;  //random walk step
double dt=0.01; //time step (=ds^2)
int mcs=500; 	//MC steps

// Variables
int N;		 //# of walkers
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
void initilize() {
	N = N0;
	for (i=1; i<N; i++)
		x[i] = rand() ; //uniform random number # in [-1,1]

	for (i=0; i<N; i++)
		V_ref += x[i]^2/2 // = <V>, this is the harmonic potential x_i^2/2
	V_ref = V_ref/N;
}

void walk() {
	Nin = N;
	Vsum = 0;
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
	}

	Vavg = Vsum/N;
	Vref = Vavg – (N-N0)/(N0*dt);
}

void data(){
	for(i=1; i<=N; i++)
		++psi[(int) (x[i] + xshift + 0.5*binsize)/binsize];
}