#include "FloatType.h"

void mp2DColl(int n, int m, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE r_omega, FLOAT_TYPE b_omega, FLOAT_TYPE control_param, FLOAT_TYPE del,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE r_A,  FLOAT_TYPE b_A, FLOAT_TYPE *r_fPert,
		FLOAT_TYPE *b_fPert, FLOAT_TYPE *weight, int *cx, int *cy);

void mp3DColl(int n, int m, int h, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *w, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE A, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl,
		FLOAT_TYPE *weight, int *cx, int *cy, int *cz, FLOAT_TYPE *f, FLOAT_TYPE r_nu, FLOAT_TYPE b_nu, FLOAT_TYPE r_alpha,
		FLOAT_TYPE b_alpha, FLOAT_TYPE *chi, FLOAT_TYPE *phi, FLOAT_TYPE *psi, FLOAT_TYPE *teta, FLOAT_TYPE *cg_w);

void createBubble(FLOAT_TYPE *x, FLOAT_TYPE *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho, int test_case);

void createBubble3D(FLOAT_TYPE *x, FLOAT_TYPE *y, FLOAT_TYPE *z, int n, int m, int h, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *phi, FLOAT_TYPE *rho, FLOAT_TYPE *f);

void updateMacroMP(int n, int m, FLOAT_TYPE *u, FLOAT_TYPE *v,FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *rho, FLOAT_TYPE control_param,
		FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE st_predicted, int test_case);

void updateMacroMP3D(int n, int m, int h, FLOAT_TYPE *u, FLOAT_TYPE *v, FLOAT_TYPE *w,FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *rho, FLOAT_TYPE control_param,
		FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE st_predicted, int *cx, int *cy, int *cz, FLOAT_TYPE *f);

void peridicBoundaries(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *u, FLOAT_TYPE *v, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *rho, int test_case);

void peridicBoundaries3D(int n, int m, int h, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho);

void streamMP(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl);

void streamMP3D(int n, int m, int h, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl, bool *stream);

void resetArrays(FLOAT_TYPE *color_gradient, int n, int m);

FLOAT_TYPE* convertArray(int n, int m, FLOAT_TYPE *arr);

void updateSurfaceTension(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE control_param,
		FLOAT_TYPE st_predicted, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, int n, int m);

FLOAT_TYPE calculateSurfaceTension(FLOAT_TYPE p_in_mean, FLOAT_TYPE p_out_mean, FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE st_predicted);

FLOAT_TYPE validateCoalescenceCase(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, int n, int m, FLOAT_TYPE radius);
