#include "Multiphase.h"
#include <math.h>
#include <stdio.h>
void mp2DColl(int *fluid, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *f, FLOAT_TYPE *fColl){

}

void createBubble(float *x, float *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho) {
	int i, j, k;
	int index, index2;
	for(i = 0; i < n; i++){
		for (j=0; j < m; j++){
			index = i*m + j;
			//printf("pow1: %f  --  pow2: %f  --  sqrt: %f\n",pow((x[index]-0.5), 2),pow((y[index]-0.5),2), sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5),2)));
			if( sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5),2)) <= radius ){
				for (k=0; k < 9; k++){
					// initialise distribution function with small, non-zero values
					r_rho[index] = r_density;
					//x + WIDTH * (y + DEPTH * z)
					index2 = k + index * 9;
					r_f[index2] = r_rho[index] * r_phi[k];
				}
				//printf("JSAHBDJSABDSHBCJHABSHDBSAHDBASHBDSHABDHSABHDBSAHDBASBDSHABDA\n");
			}
			else {
				for (k=0; k < 9; k++){
					// initialise distribution function with small, non-zero values
					b_rho[index]=b_density;
					index2 = k + index * 9;
					b_f[index2]   = b_rho[index]*b_phi[k];
				}
			}
			// initialise density
			rho[index] = r_rho[index]+b_rho[index];

			printf("%f\n", r_rho[i*m + j]);
		}
	}
}
