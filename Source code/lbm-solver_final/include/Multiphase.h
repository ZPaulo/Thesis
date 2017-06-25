#include "FloatType.h"

void mp2DColl(int *fluid,int n, int m, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE r_omega, FLOAT_TYPE b_omega, FLOAT_TYPE control_param, FLOAT_TYPE del,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE r_A,  FLOAT_TYPE b_A, FLOAT_TYPE *r_fPert, FLOAT_TYPE *b_fPert);

void createBubble(float *x, float *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho);
