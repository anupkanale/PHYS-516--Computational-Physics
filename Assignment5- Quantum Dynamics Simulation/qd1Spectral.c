/*******************************************************************************
Quantum dynamics (QD) simulation of an electron in one dimension.

USAGE

%cc -o qd1 qd1.c -lm
%qd1 < qd1.in (see qd1.h for the input-file format)
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "qd1Spectral.h"

int main(int argc, char **argv) {
	int step; /* Simulation loop iteration index */
	int i,m;

	init_param();  /* Read input parameters */
	init_prop();   /* Initialize the kinetic & potential propagators */
	init_wavefn(); /* Initialize the electron wave function */

	double refl[Nstep];
	double tran[Nstep];

	for (step=1; step<=Nstep; step++) {
		single_step(); /* Time propagation for one step, DT */

		refl[step] = 0;
		tran[step] = 0;
		for (i=0; i<NX/2; i++)
			refl[step] += (psi[2*i+0]*psi[2*i+0] + psi[2*i+1]*psi[2*i+1])*dx; 
		for (i=NX/2; i<NX; i++)
			tran[step] += (psi[2*i+0]*psi[2*i+0] + psi[2*i+1]*psi[2*i+1])*dx;
		
		// printf("%le %le %le\n", DT*step, refl[step], tran[step]);

		if (step%NECAL==0) {
			calc_energy();
			printf("%le %le %le %le\n",DT*step,ekin,epot,etot);
		}
	}

	return 0;
}

/*----------------------------------------------------------------------------*/
void init_param() {
/*------------------------------------------------------------------------------
	Initializes parameters by reading them from standard input.
------------------------------------------------------------------------------*/
	/* Read control parameters */
	scanf("%le",&LX);
	scanf("%le",&DT);
	scanf("%d",&Nstep);
	scanf("%d",&NECAL);
	scanf("%le%le%le",&X0,&S0,&E0);
	scanf("%le%le",&BH,&BW);
	scanf("%le",&EH);

	/* Calculate the mesh size */
	dx = LX/NX;
}

/*----------------------------------------------------------------------------*/
void init_prop() {
/*------------------------------------------------------------------------------
	Initializes the kinetic & potential propagators.
------------------------------------------------------------------------------*/
	int stp,s,i,up,lw,m;
	double x;

	/* Set up kinetic propagators */
	for (m=0; m<NX; m++){
		if (m<=NX/2)
			km[m] = 2*M_PI*m/LX;
		else
			km[m] = 2*M_PI*(m-NX)/LX;

		ut[2*m+0] = cos(-0.5*DT*km[m]*km[m]);
		ut[2*m+1] = sin(-0.5*DT*km[m]*km[m]);
	}

	/* Set up potential propagator */
	for (i=0; i<NX; i++) {
		x = dx*i;
		/* Construct the edge potential */
		if (i==0 || i==NX-1)
			v[i] = EH;
		/* Construct the barrier potential */
		else if (0.5*(LX-BW)<x && x<0.5*(LX+BW))
			v[i] = BH;
		else
			v[i] = 0.0;
		/* Half-step potential propagator */
		uv[2*i+0] = cos(-0.5*DT*v[i]);
		uv[2*i+1] = sin(-0.5*DT*v[i]);
	}
}

/*----------------------------------------------------------------------------*/
void init_wavefn() {
/*------------------------------------------------------------------------------
	Initializes the wave function as a traveling Gaussian wave packet.
------------------------------------------------------------------------------*/
	int sx,s;
	double x,gauss,psisq,norm_fac;

	/* Calculate the the wave function value mesh point-by-point */
	for (sx=0; sx<NX; sx++) {
		x = dx*sx-X0;
		gauss = exp(-0.25*x*x/(S0*S0));
		psi[2*sx+0] = gauss*cos(sqrt(2.0*E0)*x);
		psi[2*sx+1] = gauss*sin(sqrt(2.0*E0)*x);
	}

	/* Normalize the wave function */
	psisq=0.0;
	for (sx=0; sx<NX; sx++)
		for (s=0; s<2; s++)
			psisq += psi[2*sx+s]*psi[2*sx+s];
	psisq *= dx;
	norm_fac = 1.0/sqrt(psisq);
	for (sx=0; sx<NX; sx++)
		for (s=0; s<2; s++)
			psi[2*sx+s] *= norm_fac;
}

/*----------------------------------------------------------------------------*/
void prop( double u[]) {
	int sx;
	double wr,wi;

	for (sx=0; sx<NX; sx++) {
		wr=u[2*sx+0]*psi[2*sx+0]-u[2*sx+1]*psi[2*sx+1];
		wi=u[2*sx+0]*psi[2*sx+1]+u[2*sx+1]*psi[2*sx+0];
		psi[2*sx+0]=wr;
		psi[2*sx+1]=wi;
	}
}

/*----------------------------------------------------------------------------*/
void single_step() {
/*------------------------------------------------------------------------------
	Propagates the electron wave function for a unit time step, DT.
------------------------------------------------------------------------------*/
	int j;

	prop(uv);  /* half step potential propagation */
	four1(psi-1, (unsigned int)NX, -1);
	for (j=0; j<2*NX; j++) {
		psi[j] /= NX;
	}
	prop(ut); /* full step kinetic propagation   */
	four1(psi-1, (unsigned int)NX, 1);
	prop(uv);  /* half step potential propagation */
}

/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
void calc_energy() {
/*------------------------------------------------------------------------------
	Calculates the kinetic, potential & total energies, EKIN, EPOT & ETOT.
------------------------------------------------------------------------------*/
	int sx,s,j;
	double a,bx;

	/* Tridiagonal kinetic-energy operators */
	a =   1.0/(dx*dx);
	bx = -0.5/(dx*dx);

	four1(psi-1, (unsigned int)NX, -1);
	for (j=0; j<2*NX; j++) {
		psi[j] /= NX;
	}

	/* kinetic energy */
	ekin = 0.0;
	for (sx=0; sx<NX; sx++)
		ekin += 0.5*km[sx]*km[sx]*(psi[2*sx+0]*psi[2*sx+0] + psi[2*sx+1]*psi[2*sx+1]);
	ekin *= NX*dx;

	four1(psi-1, (unsigned int)NX, 1);
	/* Potential energy */
	epot = 0.0;
	for (sx=0; sx<NX; sx++)
		epot += v[sx]*(psi[2*sx+0]*psi[2*sx+0] + psi[2*sx+1]*psi[2*sx+1]);
	epot *= dx;

	/* Total energy */
	etot = ekin+epot;
}
