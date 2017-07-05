#include "FloatType.h"

void mp2DColl(int n, int m, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE r_omega, FLOAT_TYPE b_omega, FLOAT_TYPE control_param, FLOAT_TYPE del,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE r_A,  FLOAT_TYPE b_A, FLOAT_TYPE *r_fPert,
		FLOAT_TYPE *b_fPert, FLOAT_TYPE *weight, int *cx, int *cy);

void createBubble(FLOAT_TYPE *x, FLOAT_TYPE *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho);

void updateMacroMP(int n, int m, FLOAT_TYPE *u, FLOAT_TYPE *v,FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *rho, FLOAT_TYPE control_param,
		FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE st_predicted);

void peridicBoundaries(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE b_density, FLOAT_TYPE *u, FLOAT_TYPE *v);

void streamMP(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl);

void resetArrays(FLOAT_TYPE *color_gradient, int n, int m);

FLOAT_TYPE* convertArray(int n, int m, FLOAT_TYPE *arr);

void updateSurfaceTension(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE control_param,
		FLOAT_TYPE st_predicted, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, int n, int m);

FLOAT_TYPE calculateSurfaceTension(FLOAT_TYPE p_in_mean, FLOAT_TYPE p_out_mean, FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE st_predicted);
