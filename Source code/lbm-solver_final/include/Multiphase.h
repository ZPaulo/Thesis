#include "FloatType.h"

void mp2DColl(int *fluid, FLOAT_TYPE *rho, FLOAT_TYPE *u,
        FLOAT_TYPE *v, FLOAT_TYPE *f, FLOAT_TYPE *fColl);

void createBubble(float *x, float *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho);
