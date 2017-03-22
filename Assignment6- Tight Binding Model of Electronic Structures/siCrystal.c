#include <stdio.h>
#include <math.h>

#define NMAX 100           /* Max # of atoms */
#define NAUC 8             /* # of atoms per unit cell */
#define LCNS 5.43		   /* Lattice constant of Si (5.43 angstrom) in atomic unit */
int nAtom;                 /* # of atoms */
double r[NMAX][3];         /* r[i][0|1|2] is the x|y|z coordinate of atom i */
int InitUcell[3];          /* # of unit cells */
double RegionH[3];

double **dmatrix(int, int, int, int);
double *dvector(int, int);
void InitConf();
double SignR(double, double);
void tred2(double **, int, double *, double *);
void tqli(double *, double *, int, double **);

int main() {
	InitConf();

	double **h;	// Hamiltonian matrix
	double *d;	// Eigenvalues
	double *e;	// Work array for matrix diagonalization

	int n4;
	n4 = 4*nAtom;
	h = dmatrix(1,n4,1,n4);
	d = dvector(1,n4);
	e = dvector(1,n4);

	int ii,jj,kk;
	// Initialize Hamiltonian matrix
	for (ii=1; ii<=n4; ii++) {
		for (jj=1; jj<=n4; jj++) {
			h[ii][jj] = 0;
		}
	}

	/* Populate the Hamiltonian matrix */
	double r0 = 2.360352, n = 2, Es = -5.25, Ep = 1.2; // Constants

	// index-- ss\sigma sp\sigma pp\sigma pp\pi
	double hLambda_r0[] = {-2.038, 1.745, 2.75, -1.075};
	double nLambda[] = {9.5, 8.5, 7.5, 7.5};
	double rLambda[] = {3.4, 3.55, 3.7, 3.7};

	double hLambda[] = {0,0,0,0};
	double rx,ry,rz,rMag,dx,dy,dz;
	double rij[] = {0,0,0};
	double expterm;

	for (jj=0; jj<nAtom; jj++) {
		for (ii=0; ii<nAtom; ii++) {
		  if (ii==jj) {
			// Diagonal elements
			h[1+4*jj][1+4*ii] = Es;
			h[2+4*jj][2+4*ii] = Ep;
			h[3+4*jj][3+4*ii] = Ep;
			h[4+4*jj][4+4*ii] = Ep;

			// Off-diagonal elements are zero
		  }

		  else {
			// Calc rij with min image convention
			for (kk=0; kk<3; kk++) {
				rij[kk] = r[ii][kk] - r[jj][kk];
				/* Chooses the nearest image */
				rij[kk] = rij[kk] - SignR(RegionH[kk],rij[kk]-RegionH[kk]) - SignR(RegionH[kk],rij[kk]+RegionH[kk]);
			}

			rx = rij[0];
			ry = rij[1];
			rz = rij[2];
			rMag = sqrt(rx*rx + ry*ry + rz*rz); 

			dx = rx/rMag; dy = ry/rMag; dz = rz/rMag; // Unit vector

			expterm = n*(-pow((rMag/rLambda[0]), nLambda[0]) + pow((r0/rLambda[0]),nLambda[0]));
			hLambda[0] = hLambda_r0[0]* pow((r0/rMag),n)* exp(expterm); // h_ss\sigma

			expterm = n*(-pow((rMag/rLambda[1]), nLambda[1]) + pow((r0/rLambda[1]),nLambda[1]));
			hLambda[1] = hLambda_r0[1]* pow((r0/rMag),n)* exp(expterm); // h_sp\sigma

			expterm = n*(-pow((rMag/rLambda[2]), nLambda[2]) + pow((r0/rLambda[2]),nLambda[2]));
			hLambda[2] = hLambda_r0[2]* pow((r0/rMag),n)* exp(expterm); // h_ss\sigma

			expterm = n*(-pow((rMag/rLambda[3]), nLambda[3]) + pow((r0/rLambda[3]),nLambda[3]));
			hLambda[3] = hLambda_r0[3]* pow((r0/rMag),n)* exp(expterm); // h_pp\pi

			// Diagonal elements
			h[1+4*jj][1+4*ii] = hLambda[0];
			h[2+4*jj][2+4*ii] = dx*dx*hLambda[2] + (1-dx*dx)*hLambda[3];
			h[3+4*jj][3+4*ii] = dy*dy*hLambda[2] + (1-dy*dy)*hLambda[3];
			h[4+4*jj][4+4*ii] = dz*dz*hLambda[2] + (1-dz*dz)*hLambda[3];

			// Populate off-Diagonal elements
			h[1+4*jj][2+4*ii] = dx*hLambda[1]; h[2+4*jj][1+4*ii] = -h[1+4*jj][2+4*ii]; 
			h[1+4*jj][3+4*ii] = dy*hLambda[1]; h[3+4*jj][1+4*ii] = -h[1+4*jj][3+4*ii];
			h[1+4*jj][4+4*ii] = dz*hLambda[1]; h[4+4*jj][1+4*ii] = -h[1+4*jj][4+4*ii];

			h[2+4*jj][3+4*ii] = dx*dy*(hLambda[2] - hLambda[3]); h[3+4*jj][2+4*ii] = h[2+4*jj][3+4*ii];
			h[2+4*jj][4+4*ii] = dx*dz*(hLambda[2] - hLambda[3]); h[4+4*jj][2+4*ii] = h[2+4*jj][4+4*ii];

			h[3+4*jj][4+4*ii] = dy*dz*(hLambda[2] - hLambda[3]); h[4+4*jj][3+4*ii] = h[3+4*jj][4+4*ii];
		  }
		}
	}

	/* Diagonalize the Hamiltonian matrix */
	tred2(h,n4,d,e);
	tqli(d,e,n4,h);

	// Print out the Hamiltonian matrix
	for (ii=1; ii<=n4; ii++) {
		printf("\n");
		for (jj=1; jj<=n4; jj++) {
			printf("%12le ", h[ii][jj]);
		}
	}
}


double SignR(double v,double x) {
	if (x > 0) return v;
	else return -v;
	}

void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to diamond lattice positions.  
------------------------------------------------------------------------------*/
	double gap[3];      /* Unit cell size */
	double c[3];
	int j,k,nX,nY,nZ;
	/* Atom positions in a unit diamond crystalline unit cell */
	double origAtom[NAUC][3] = {{0.0, 0.0, 0.0 }, {0.0, 0.5, 0.5 },
                                {0.5, 0.0, 0.5 }, {0.5, 0.5, 0.0 },
                                {0.25,0.25,0.25}, {0.25,0.75,0.75},
                                {0.75,0.25,0.75}, {0.75,0.75,0.25}};

	/* Read the # of unit cells in the x, y & z directions */
	scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
	
	for (k=0; k<3; k++){
		RegionH[k] = 0.5*LCNS*InitUcell[k];
	}

	/* Sets up a diamond lattice */
	for (k=0; k<3; k++) gap[k] = LCNS;
	nAtom = 0;
	for (nZ=0; nZ<InitUcell[2]; nZ++) {
		c[2] = nZ*gap[2];
		for (nY=0; nY<InitUcell[1]; nY++) {
			c[1] = nY*gap[1];
			for (nX=0; nX<InitUcell[0]; nX++) {
				c[0] = nX*gap[0];
				for (j=0; j<NAUC; j++) {
					for (k=0; k<3; k++)
						r[nAtom][k] = c[k] + gap[k]*origAtom[j][k];
					++nAtom;
				}
			}
		}
	}
}
