#include "Multiphase.h"
#include <math.h>
#include <stdio.h>
void mp2DColl(int *fluid,int n, int m, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE r_omega, FLOAT_TYPE b_omega, FLOAT_TYPE control_param, FLOAT_TYPE del,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE r_A,  FLOAT_TYPE b_A, FLOAT_TYPE *r_fPert, FLOAT_TYPE *b_fPert){

	FLOAT_TYPE cu1, cu2;

	int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	FLOAT_TYPE weight[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
	FLOAT_TYPE cosin[9] = {0.0};
	FLOAT_TYPE chi;
	FLOAT_TYPE r_omega_temp, b_omega_temp;
	FLOAT_TYPE a1      =   2.0 * r_omega * b_omega/(r_omega+b_omega);
	FLOAT_TYPE a2      =   2.0 * (r_omega - a1) / del;
	FLOAT_TYPE a3      =   -a2 / (2.0 * del);
	FLOAT_TYPE a4      =   2.0 * (a1 - b_omega) / del;
	FLOAT_TYPE a5      =   a4 / (2.0 * del);
	FLOAT_TYPE color_gradient_norm;
	FLOAT_TYPE k_r, k_b, k_k;
	FLOAT_TYPE norm_c;
	FLOAT_TYPE prod_c_g;
	FLOAT_TYPE r_pert, b_pert;
	FLOAT_TYPE r_feq, b_feq;
	FLOAT_TYPE fn05;
	int index, index9, temp_index;
	for (int i=0;i < n; i++){
		for (int j=0; j < m; j++){
			// temporary variable 1
			index = i*m + j;
			cu1 = u[index]*u[index] + v[index]*v[index];

			for (int k=0; k<9; k++){
				// temporary variable 2
				cu2 = u[index]*cx[k] + v[index]*cy[k];

				index9 = k + index * 9;
				// calculate equilibrium distribution function


				// calculate color gradient - 4th order
				if (k!=0){ // the rest node (k=0) does not contribute to the color gradient

					if (i!=0 && j!=0 && i!=(n-1) && j!=(m-1)){ // Interior points - In the boundary it is calculated by "mirroring" the density
						temp_index = (i + cx[k]) * m + j + cy[k];
						color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k];
						color_gradient[index * 2 + 1] += (r_rho[temp_index] - b_rho[temp_index]) * cy[k];
					}
					else if (j==(m-1) && i!=0 && i!=(n-1)) {// north boundary
						temp_index = (i + cx[k]) * m + j - abs(cy[k]);
						color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k];
						color_gradient[index * 2 + 1] = 0;
					}
					else if (j==0 && i!=0 && i!=(n-1)){  // south boundary
						temp_index = (i + cx[k]) * m + j + abs(cy[k]);
						color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k];
						color_gradient[index * 2 + 1] = 0;
					}
					else if (i==(n-1) && j!=0 && j!=(m-1)){  // east boundary
						temp_index = (i - abs(cx[k])) * m + j + cy[k];
						color_gradient[index * 2] = 0;
						color_gradient[index * 2 + 1] +=  (r_rho[temp_index] - b_rho[temp_index]) * cy[k];
					}
					else if (i==0 && j!=0 && j!=(m-1)){ //  west boundary
						temp_index = (i + abs(cx[k])) * m + j + cy[k];
						color_gradient[index * 2] = 0;
						color_gradient[index * 2 + 1] +=  (r_rho[temp_index] - b_rho[temp_index]) * cy[k];
					}
				}
			}

			// relaxation parameter to choose a proper omega at the interface
			if (r_omega != b_omega){
				chi=(r_rho[index] - b_rho[index])/rho[index];
				if(chi >= -control_param && chi <= control_param){
					if (chi > del)
						r_omega_temp=r_omega;
					else if (chi <= del && chi > 0)
						r_omega_temp=a1 + a2 * chi + a3 * chi * chi;
					else if (chi <= 0 && chi >= -del)
						r_omega_temp=a1 + a4 * chi + a5 * chi * chi;
					else if (chi < -del)
						r_omega_temp=b_omega;
				}
			}
			else
				r_omega_temp=r_omega;

			b_omega_temp=r_omega_temp;

			// invariable quantities
			color_gradient_norm = sqrt(pow(color_gradient[index * 2],2) + pow(color_gradient[index * 2 + 1],2));
			k_r=r_rho[index]/rho[index];
			k_b=b_rho[index]/rho[index];
			k_k= beta * r_rho[index] * b_rho[index]/(pow(rho[index],2));

			for (int k=0;k<9;k++){
				if (color_gradient_norm > g_limit){
					if (k!=0){
						norm_c= sqrt(pow(cx[k],2)+pow(cy[k],2));
						cosin[k]=(cx[k]*color_gradient[index * 2]+cy[k]*color_gradient[index * 2 + 1]) / (color_gradient_norm*norm_c);
					}
					else
						cosin[k]=0.0;
					// calculate perturbation terms
					prod_c_g=cx[k]*color_gradient[index * 2]+cy[k]*color_gradient[index * 2 + 1];
					r_pert=0.5*r_A*color_gradient_norm*(weight[k]*pow((prod_c_g/color_gradient_norm),2)-w_pert[k]);
					b_pert=0.5*b_A*color_gradient_norm*(weight[k]*pow((prod_c_g/color_gradient_norm),2)-w_pert[k]);

				}
				else{
					// ther perturbation terms are null
					r_pert=0;
					b_pert=0;
				}
				// calculate updated distribution function
				index9 = k + index * 9;
				r_feq = r_rho[index] * (r_phi[k] + weight[k] * (3 * cu2 + 4.5 * cu2 * cu2 - 1.5 * cu1 * cu1));
				b_feq = b_rho[index] * (b_phi[k] + weight[k] * (3 * cu2 + 4.5 * cu2 * cu2 - 1.5 * cu1 * cu1));

				r_fPert[index9] = r_omega_temp*r_feq + (1-r_omega_temp)*r_f[index9]+r_pert;
				b_fPert[index9] = b_omega_temp*b_feq + (1-b_omega_temp)*b_f[index9]+b_pert;
				fn05 = r_fPert[index9] + b_fPert[index9];
				// perform recolor step
				r_fPert[index9]=k_r*fn05+k_k*cosin[k]*(r_rho[index]*r_phi[k]+b_rho[index]*b_phi[k]);
				b_fPert[index9]=k_b*fn05-k_k*cosin[k]*(r_rho[index]*r_phi[k]+b_rho[index]*b_phi[k]);
			}
		}
	}



}

void createBubble(float *x, float *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho) {
	int i, j, k;
	int index, index2;
	for(i = 0; i < n; i++){
		for (j=0; j < m; j++){
			index = i*m + j;
			if( sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5),2)) <= radius ){
				for (k=0; k < 9; k++){
					// initialise distribution function with small, non-zero values
					r_rho[index] = r_density;
					//x + WIDTH * (y + DEPTH * z)
					index2 = k + index * 9;
					r_f[index2] = r_rho[index] * r_phi[k];
				}
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
		}
	}
}
