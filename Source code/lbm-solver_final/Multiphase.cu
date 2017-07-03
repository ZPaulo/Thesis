#include "Multiphase.h"
#include <math.h>
#include <stdio.h>
#include "ArrayUtils.h"

void mp2DColl(int n, int m, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE r_omega, FLOAT_TYPE b_omega, FLOAT_TYPE control_param, FLOAT_TYPE del,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE r_A,  FLOAT_TYPE b_A, FLOAT_TYPE *r_fPert, FLOAT_TYPE *b_fPert,
		FLOAT_TYPE *weight, int *cx, int *cy){

	FLOAT_TYPE cu1, cu2, r_CollPert, b_CollPert;
	FLOAT_TYPE cosin;
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
	for (int j=0; j < m; j++){
		for (int i=0;i < n; i++){
			index = j*n + i;


			color_gradient[index * 2] = 0;
			color_gradient[index * 2 + 1] = 0;
			for (int k=0; k<9; k++){

				// calculate color gradient - 4th order

				if (i!=0 && j!=0 && i!=(n-1) && j!=(m-1)){ // Interior points - In the boundary it is calculated by "mirroring" the density
					temp_index = (j + cy[k]) * n + i + cx[k];
					color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k];
					color_gradient[index * 2 + 1] += (r_rho[temp_index] - b_rho[temp_index]) * cy[k];
				}
				else if (j==(m-1) && i!=0 && i!=(n-1)) {// north boundary
					temp_index = (j - abs(cy[k])) * n + i + cx[k];
					color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k];
					color_gradient[index * 2 + 1] = 0;
				}
				else if (j==0 && i!=0 && i!=(n-1)){  // south boundary
					temp_index = (j + abs(cy[k])) * n + i + cx[k];
					color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k];
					color_gradient[index * 2 + 1] = 0;
				}
				else if (i==(n-1) && j!=0 && j!=(m-1)){  // east boundary
					temp_index = (j + cy[k]) * n + i - abs(cx[k]);
					color_gradient[index * 2] = 0;
					color_gradient[index * 2 + 1] +=  (r_rho[temp_index] - b_rho[temp_index]) * cy[k];
				}
				else if (i==0 && j!=0 && j!=(m-1)){ //  west boundary
					temp_index = (j + cy[k]) * n + i + abs(cx[k]);
					color_gradient[index * 2] = 0;
					color_gradient[index * 2 + 1] +=  (r_rho[temp_index] - b_rho[temp_index]) * cy[k];
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

			cu1 = u[index]*u[index] + v[index]*v[index];

			// invariable quantities
			color_gradient_norm = sqrt(pow(color_gradient[index * 2],2) + pow(color_gradient[index * 2 + 1],2));
			k_r=r_rho[index]/rho[index];
			k_b=b_rho[index]/rho[index];
			k_k= beta * r_rho[index] * b_rho[index]/(pow(rho[index],2));
			for (int k=0;k<9;k++){
				if (color_gradient_norm > g_limit){
					prod_c_g=cx[k]*color_gradient[index * 2]+cy[k]*color_gradient[index * 2 + 1];
					if (k!=0){
						norm_c= sqrt(pow(cx[k],2)+pow(cy[k],2));
						cosin= prod_c_g / (color_gradient_norm*norm_c);
					}
					else
						cosin=0.0;
					// calculate perturbation terms

					r_pert=0.5*r_A*color_gradient_norm*(weight[k]*pow((prod_c_g/color_gradient_norm),2)-w_pert[k]);
					b_pert=0.5*b_A*color_gradient_norm*(weight[k]*pow((prod_c_g/color_gradient_norm),2)-w_pert[k]);

				}
				else{
					// ther perturbation terms are null
					r_pert=0.0;
					b_pert=0.0;
				}


				cu2 = u[index]*cx[k] + v[index]*cy[k];
				// calculate equilibrium distribution function
				r_feq = r_rho[index] * (r_phi[k] + weight[k] * (3 * cu2 + 4.5 * cu2 * cu2 - 1.5 * cu1));
				b_feq = b_rho[index] * (b_phi[k] + weight[k] * (3 * cu2 + 4.5 * cu2 * cu2 - 1.5 * cu1));

				index9 = i + j * n + k * m * n;
				// calculate updated distribution function
				r_CollPert = r_omega_temp*r_feq + (1-r_omega_temp)*r_f[index9]+r_pert;
				b_CollPert = b_omega_temp*b_feq + (1-b_omega_temp)*b_f[index9]+b_pert;

				fn05 = r_CollPert + b_CollPert;
				//				// perform recolor step
				r_fPert[index9]=k_r*fn05+k_k*cosin*(r_rho[index]*r_phi[k]+b_rho[index]*b_phi[k]);
				b_fPert[index9]=k_b*fn05-k_k*cosin*(r_rho[index]*r_phi[k]+b_rho[index]*b_phi[k]);
			}
		}
	}



}

void createBubble(FLOAT_TYPE *x, FLOAT_TYPE *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho) {
	int i, j, k;
	int index, index2;
	for (j=0; j < m; j++){
		for(i = 0; i < n; i++){
			index = j * n + i;

			if( sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5),2)) <= radius ){
				r_rho[index] = r_density;
				for (k=0; k < 9; k++){
					// initialise distribution function with small, non-zero values
					index2 = i + j * n + k * m * n;
					r_f[index2] = r_rho[index] * r_phi[k];
				}
			}
			else {
				b_rho[index]=b_density;
				for (k=0; k < 9; k++){
					// initialise distribution function with small, non-zero values
					index2 = i + j * n + k * m * n;
					b_f[index2]   = b_rho[index]*b_phi[k];
				}
			}
			// initialise density
			rho[index] = r_rho[index]+b_rho[index];
		}
	}
}

void updateMacroMP(int n, int m, FLOAT_TYPE *u, FLOAT_TYPE *v,FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *rho, FLOAT_TYPE control_param,
		FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE st_predicted){

	int index_aux1=0;
	int index_aux2=0;
	FLOAT_TYPE p_in=0.0;
	FLOAT_TYPE p_out=0.0;
	FLOAT_TYPE u_cum, v_cum;
	FLOAT_TYPE r_sum, b_sum;
	FLOAT_TYPE chi;
	int index, index9;
	int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	FLOAT_TYPE st_laplace;
	// Density and Velocity
	for (int j=1; j < m - 1;j++){
		for (int i=1; i < n - 1; i++){
			// auxiliar variables
			u_cum=0.0;
			v_cum=0.0;

			// densities
			index = j * n + i;
			r_sum = 0.0;
			b_sum = 0.0;
			for(int k = 0; k < 9; k++){
				index9 = i + j * n + k * m * n;
				r_sum += r_f[index9];
				b_sum += b_f[index9];
			}
			r_rho[index] = r_sum;
			b_rho[index]= b_sum;
			rho[index] = r_rho[index]+b_rho[index];

			// p_in and p_out for the surface tension
			chi=(r_rho[index]-b_rho[index])/rho[index];
//			printf("chi "FLOAT_FORMAT" ",chi);
//			printf("control "FLOAT_FORMAT" \n", control_param);
			if (chi >= control_param){
				index_aux1++;
				p_in += r_rho[index];
			}
			else if (chi < control_param){
				index_aux2++;
				p_out+=b_rho[index];
			}

			// velocities
			for (int k=0; k < 9; k++){
				index9 = i + j * n + k * m * n;
				u_cum += (r_f[index9]+b_f[index9])*cx[k];
				v_cum += (r_f[index9]+b_f[index9])*cy[k];
			}
			u[index]   = u_cum/rho[index];
			v[index]  = v_cum/rho[index];

		}
	}

	// Calculate surface tension
	//printf("%f vs %f\n", p_in, p_out);
	p_in=(3.0/5.0)*(1.0-r_alpha)*p_in/index_aux1;      // pressure average inside the bubble
	p_out=(3.0/5.0)*(1.0-b_alpha)*p_out/index_aux2;   // pressure average outside the bubble
	st_laplace=bubble_radius*(p_in-p_out);

	st_error[iteration]=abs(st_predicted-st_laplace)/(st_predicted)*100.0;
}

void peridicBoundaries(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE b_density, FLOAT_TYPE *u, FLOAT_TYPE *v){

	int index_end, index_start;
	int jn = m-1;
	int js = 0;
	int ie = n-1;
	int iw = 0;

	for (int i=1; i < n-1; i++){
		// north boundary
		index_end = jn * n + i;
		index_start = js * n + i;

		r_f[index_end + 4 * m * n] = r_f[index_start + 4 * m * n];
		r_f[index_end + 7 * m * n] = r_f[index_start + 7 * m * n];
		r_f[index_end + 8 * m * n] = r_f[index_start + 8 * m * n];

		b_f[index_end + 4 * m * n] = b_f[index_start + 4 * m * n];
		b_f[index_end + 7 * m * n] = b_f[index_start + 7 * m * n];
		b_f[index_end + 8 * m * n] = b_f[index_start + 8 * m * n];

		// macroscopic boundary conditions
//		r_rho[index_end] = 0;
//		b_rho[index_end] = b_density;
//		u[index_end]   = 0;
//		v[index_end]   = 0;

		//south boundary
		r_f[index_start + 2 * m * n] = r_f[index_end + 2 * m * n];
		r_f[index_start + 5 * m * n] = r_f[index_end + 5 * m * n];
		r_f[index_start + 6 * m * n] = r_f[index_end + 6 * m * n];

		b_f[index_start + 2 * m * n] = b_f[index_end + 2 * m * n];
		b_f[index_start + 5 * m * n] = b_f[index_end + 5 * m * n];
		b_f[index_start + 6 * m * n] = b_f[index_end + 6 * m * n];

//		r_rho[index_start] = 0;
//		b_rho[index_start] = b_density;
//		u[index_start]   = 0;
//		v[index_start]   = 0;
	}



	for (int j=1; j < m-1; j++){
		// east boundary
		index_end = j*n + ie;
		index_start = j*n + iw;

		r_f[index_end + 3 * m * n] = r_f[index_start + 3 * m * n];
		r_f[index_end + 7 * m * n] = r_f[index_start + 7 * m * n];
		r_f[index_end + 6 * m * n] = r_f[index_start + 6 * m * n];

		b_f[index_end + 3 * m * n] = b_f[index_start + 3 * m * n];
		b_f[index_end + 7 * m * n] = b_f[index_start + 7 * m * n];
		b_f[index_end + 6 * m * n] = b_f[index_start + 6 * m * n];

		//macroscopic boundary conditions
//		r_rho[index_end] = 0;
//		b_rho[index_end] = b_density;
//		u[index_end]   = 0;
//		v[index_end]   = 0;

		// west boundary
		r_f[index_start + 1 * m * n] = r_f[index_end + 1 * m * n];
		r_f[index_start + 5 * m * n] = r_f[index_end + 5 * m * n];
		r_f[index_start + 8 * m * n] = r_f[index_end + 8 * m * n];

		b_f[index_start + 1 * m * n] = b_f[index_end + 1 * m * n];
		b_f[index_start + 5 * m * n] = b_f[index_end + 5 * m * n];
		b_f[index_start + 8 * m * n] = b_f[index_end + 8 * m * n];

//		r_rho[index_start] = 0;
//		b_rho[index_start] = b_density;
//		u[index_start]   = 0;
//		v[index_start]   = 0;

	}

	// north-east corner
	r_f[(jn*n+ie) + 3 * m * n] = r_f[(jn*n+iw) + 3 * m * n];
	r_f[(jn*n+ie) + 4 * m * n] = r_f[(js*n+ie) + 4 * m * n];
	r_f[(jn*n+ie) + 7 * m * n] = r_f[(js*n+iw) + 7 * m * n];

	b_f[(jn*n+ie) + 3 * m * n] = b_f[(jn*n+iw) + 3 * m * n];
	b_f[(jn*n+ie) + 4 * m * n] = b_f[(js*n+ie) + 4 * m * n];
	b_f[(jn*n+ie) + 7 * m * n] = b_f[(js*n+iw) + 7 * m * n];

	//	FLOAT_TYPE sum_r = 0.0;
	//	FLOAT_TYPE sum_b = 0.0;
	//	for(int i = 0; i < 9; i++){
	//		sum_r += r_f[(jn*n+ie)*9 + i];
	//		sum_b += b_f[(jn*n+ie)*9 + i];
	//	}
	//
	//	r_rho[jn*n+ie] = sum_r;
	//	b_rho[jn*n+ie] = sum_b;

//	r_rho[jn*n+ie] = 0;
//	b_rho[jn*n+ie] = b_density;
//
//	u[jn*n+ie]   = 0;
//	v[jn*n+ie]   = 0;

	// north-west corner
	r_f[(jn*n+iw) + 1 * m * n] = r_f[(jn*n+ie) + 1 * m * n];
	r_f[(jn*n+iw) + 4 * m * n] = r_f[(js*n+iw) + 4 * m * n];
	r_f[(jn*n+iw) + 8 * m * n] = r_f[(js*n+ie) + 8 * m * n];

	b_f[(jn*n+iw) + 1 * m * n] = b_f[(jn*n+ie) + 1 * m * n];
	b_f[(jn*n+iw) + 4 * m * n] = b_f[(js*n+iw) + 4 * m * n];
	b_f[(jn*n+iw) + 8 * m * n] = b_f[(js*n+ie) + 8 * m * n];

	//	sum_r = 0.0;
	//	sum_b = 0.0;
	//	for(int i = 0; i < 9; i++){
	//		sum_r += r_f[(jn*n+iw)*9 + i];
	//		sum_b += b_f[(jn*n+iw)*9 + i];
	//	}
	//
	//	r_rho[jn*n+iw] = sum_r;
	//	b_rho[jn*n+iw] = sum_b;

//	r_rho[jn*n+iw] = 0;
//	b_rho[jn*n+iw] = b_density;
//
//	u[jn*n+iw]   = 0;
//	v[jn*n+iw]   = 0;

	// south-east corner
	r_f[(js*n+ie) + 2 * m * n] = r_f[(jn*n+ie) + 2 * m * n];
	r_f[(js*n+ie) + 3 * m * n] = r_f[(js*n+iw) + 3 * m * n];
	r_f[(js*n+ie) + 6 * m * n] = r_f[(jn*n+iw) + 6 * m * n];

	b_f[(js*n+ie) + 2 * m * n] = b_f[(jn*n+ie) + 2 * m * n];
	b_f[(js*n+ie) + 3 * m * n] = b_f[(js*n+iw) + 3 * m * n];
	b_f[(js*n+ie) + 6 * m * n] = b_f[(jn*n+iw) + 6 * m * n];

	//	sum_r = 0.0;
	//	sum_b = 0.0;
	//	for(int i = 0; i < 9; i++){
	//		sum_r += r_f[(js*n+ie)*9 + i];
	//		sum_b += b_f[(js*n+ie)*9 + i];
	//	}
	//
	//	r_rho[js*n+ie] = sum_r;
	//	b_rho[js*n+ie] = sum_b;

//	r_rho[js*n+ie] = 0;
//	b_rho[js*n+ie] = b_density;
//
//	u[js*n+ie]   = 0;
//	v[js*n+ie]   = 0;


	// south-west corner
	r_f[(js*n+iw) + 2 * m * n] = r_f[(jn*n+iw) + 2 * m * n];
	r_f[(js*n+iw) + 1 * m * n] = r_f[(js*n+ie) + 1 * m * n];
	r_f[(js*n+iw) + 5 * m * n] = r_f[(jn*n+ie) + 5 * m * n];

	b_f[(js*n+iw) + 2 * m * n] = b_f[(jn*n+iw) + 2 * m * n];
	b_f[(js*n+iw) + 1 * m * n] = b_f[(js*n+ie) + 1 * m * n];
	b_f[(js*n+iw) + 5 * m * n] = b_f[(jn*n+ie) + 5 * m * n];

	//	sum_r = 0.0;
	//	sum_b = 0.0;
	//	for(int i = 0; i < 9; i++){
	//		sum_r += r_f[(js*n+iw)*9 + i];
	//		sum_b += b_f[(js*n+iw)*9 + i];
	//	}
	//
	//	r_rho[js*n+iw] = sum_r;
	//	b_rho[js*n+iw] = sum_b;

//	r_rho[js*n+iw] = 0;
//	b_rho[js*n+iw] = b_density;
//
//	u[js*n+iw]   = 0;
//	v[js*n+iw]  = 0;
}

void streamMP(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl){
	// stream on interior first
	int index,i,j;
	for (j=1;j < m-1;j++){
		for (i=1; i < n-1; i++){
			index = j*n+i;
			r_f[index] = r_fColl[index];
			r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
			r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
			r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
			r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
			r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];
			r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];
			r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];
			r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

			b_f[index] = b_fColl[index];
			b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
			b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
			b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
			b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
			b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];
			b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];
			b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];
			b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];
		}
	}
	for (i=1; i < n-1; i++){
		//north boundary
		j = m-1;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
		r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
		r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
		r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];
		r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
		b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
		b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
		b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];
		b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];

		//South boundary
		j = 0;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
		r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
		r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
		r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];
		r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
		b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
		b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
		b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];
		b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];
	}

	for (j=1;j < m-1;j++){
		//east
		i = n-1;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
		r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
		r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
		r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];
		r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
		b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
		b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
		b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];
		b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];

		//west
		i = 0;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
		r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
		r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
		r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];
		r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
		b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
		b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
		b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];
		b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];
	}

	// north-east corner
	i=n-1; j=m-1;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
	r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
	r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
	b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
	b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];

	//north-west corner
	i=0; j=m-1;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
	r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
	r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
	b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
	b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];

	// south-east corner
	i=n-1; j=0;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
	r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
	r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
	b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
	b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];

	// south-west corner
	i=0; j=0;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
	r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
	r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
	b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
	b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];

}

void resetArrays(FLOAT_TYPE *color_gradient, int n, int m){
	for(int i = 0; i < m * n *2; i++){
		color_gradient[i] = 0.0;
	}
}

FLOAT_TYPE* convertArray(int n, int m, FLOAT_TYPE *arr){
	FLOAT_TYPE *result = createHostArrayFlt(n*m, ARRAY_NONE);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			result[j*n+i] = arr[i*m+j];
		}
	}

	return result;
}
