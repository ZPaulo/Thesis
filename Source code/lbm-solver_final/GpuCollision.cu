/**
 * Collision model
 * @file GpuCollision.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */

#include "GpuFunctions.h"
#include "GpuConstants.h"

/**
 * @brief Compute the equilibrum distribution without the collision frequency
 * @details ...
 * @todo doc: insert name or reference, explain method
 *
 * @param u,v velocity
 * @param uc,vc velocity component, see: #cx_d #cy_d
 * @param rho density
 * @param w lattice weight, see: #w_d
 * @return equlibrum distribution function
 */
__device__ FLOAT_TYPE feqc2D(FLOAT_TYPE u, FLOAT_TYPE uc, FLOAT_TYPE v, FLOAT_TYPE vc, FLOAT_TYPE rho, FLOAT_TYPE w)
{
	FLOAT_TYPE vec = u*uc + v*vc;
	return rho * w * (1. + 3.*vec + 4.5*vec*vec - 1.5*(u*u + v*v));
}

__device__ FLOAT_TYPE feqc2DCG(FLOAT_TYPE u, FLOAT_TYPE uc, FLOAT_TYPE v, FLOAT_TYPE vc, FLOAT_TYPE rho, FLOAT_TYPE w, FLOAT_TYPE phi)
{
	FLOAT_TYPE vec = u*uc + v*vc;
	return rho * (phi + w * (3.*vec + 4.5*vec*vec - 1.5*(u*u + v*v)));
}

__device__ FLOAT_TYPE feqc3D(FLOAT_TYPE u, FLOAT_TYPE uc, FLOAT_TYPE v, FLOAT_TYPE vc, FLOAT_TYPE w, FLOAT_TYPE wc, FLOAT_TYPE rho, FLOAT_TYPE weight)
{
	FLOAT_TYPE vec = u*uc + v*vc + w*wc;

	return rho * weight * (1. + 3.*vec + 4.5*vec*vec - 1.5*(u*u + v*v+ w*w));
}

__device__ FLOAT_TYPE calculateColorGradientY(FLOAT_TYPE *r_rho_d, FLOAT_TYPE *b_rho_d,
		int cg_dir_d, int index){

	FLOAT_TYPE res = 0.0;
	int n = length_d;

	//int cx2D[9] = { 0, 1, 0, -1,  0, 1, -1, -1,  1 };
	//int cy2D[9] = { 0, 0, 1,  0, -1, 1,  1, -1, -1 };
	FLOAT_TYPE cg_w1 = 4.0/12.0;
	FLOAT_TYPE cg_w2 = 1.0/12.0;
	switch (cg_dir_d) {
	case 0:
		res += (r_rho_d[index + n] - b_rho_d[index + n])  *  cg_w1;
		res += (r_rho_d[index - n] - b_rho_d[index - n])  * -cg_w1;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1])  *  cg_w2;
		res += (r_rho_d[index + n - 1] - b_rho_d[index + n - 1]) * cg_w2;
		res += (r_rho_d[index - n - 1] - b_rho_d[index - n - 1])   * -cg_w2;
		res += (r_rho_d[index - n + 1] - b_rho_d[index - n + 1])  *  -cg_w2;
		break;
	case 3: //EAST
		res += (r_rho_d[index + n] - b_rho_d[index + n])  *  cg_w1;
		res += (r_rho_d[index - n] - b_rho_d[index - n])  * -cg_w1;
		res += (r_rho_d[index + n - 1] - b_rho_d[index + n - 1])  *  cg_w2;
		res += (r_rho_d[index + n - 1] - b_rho_d[index + n - 1]) * cg_w2;
		res += (r_rho_d[index - n - 1] - b_rho_d[index - n - 1])   * -cg_w2;
		res += (r_rho_d[index - n - 1] - b_rho_d[index - n - 1])  *  -cg_w2;
		break;
	case 4: //WEST
		res += (r_rho_d[index + 1] - b_rho_d[index + 1])  *  cg_w1;
		res += (r_rho_d[index - 1] - b_rho_d[index - 1])  * -cg_w1;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1])  *  cg_w2;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1]) * -cg_w2;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1])   * -cg_w2;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1])  *  cg_w2;
		break;
	default:
		break;
	}

	return res;
}

__device__ FLOAT_TYPE calculateColorGradientX(FLOAT_TYPE *r_rho_d, FLOAT_TYPE *b_rho_d,
		int cg_dir_d, int index){

	FLOAT_TYPE res = 0.0;
	int n = length_d;

	//int cx2D[9] = { 0, 1, 0, -1,  0, 1, -1, -1,  1 };
	//int cy2D[9] = { 0, 0, 1,  0, -1, 1,  1, -1, -1 };
	//cg weights = {0, 4, 4, 4, 4, 1, 1, 1, 1}/12
	FLOAT_TYPE cg_w1 = 4.0/12.0;
	FLOAT_TYPE cg_w2 = 1.0/12.0;
	switch (cg_dir_d) {
	case 0:
		res += (r_rho_d[index + 1] - b_rho_d[index + 1])  * cg_w1;
		res += (r_rho_d[index - 1] - b_rho_d[index - 1])  * -cg_w1;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1])  *  cg_w2;
		res += (r_rho_d[index + n - 1] - b_rho_d[index + n - 1]) * -cg_w2;
		res += (r_rho_d[index - n - 1] - b_rho_d[index - n - 1])   * -cg_w2;
		res += (r_rho_d[index - n + 1] - b_rho_d[index - n + 1])  *  cg_w2;
		break;
	case 1: //NORTH
		res += (r_rho_d[index + 1] - b_rho_d[index + 1])  *  cg_w1;
		res += (r_rho_d[index - 1] - b_rho_d[index - 1])  * -cg_w1;
		res += (r_rho_d[index - n + 1] - b_rho_d[index - n + 1])  *  cg_w2;
		res += (r_rho_d[index - n - 1] - b_rho_d[index - n - 1]) * -cg_w2;
		res += (r_rho_d[index - n - 1] - b_rho_d[index - n - 1])   * -cg_w2;
		res += (r_rho_d[index - n + 1] - b_rho_d[index - n + 1])  *  cg_w2;
		break;
	case 2: //SOUTH
		res += (r_rho_d[index + 1] - b_rho_d[index + 1])  *  cg_w1;
		res += (r_rho_d[index - 1] - b_rho_d[index - 1])  * -cg_w1;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1])  *  cg_w2;
		res += (r_rho_d[index + n - 1] - b_rho_d[index + n - 1]) * -cg_w2;
		res += (r_rho_d[index + n - 1] - b_rho_d[index + n - 1])   * -cg_w2;
		res += (r_rho_d[index + n + 1] - b_rho_d[index + n + 1])  *  cg_w2;
		break;
	default:
		break;
	}

	return res;
}

__global__ void gpuCollBgkw2D(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d,
		FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;

	int ms = depth_d*length_d;
	FLOAT_TYPE r, u, v;
	if (ind < ms && fluid_d[ind] == 1)
	{
		u =   u_d[ind];
		v =   v_d[ind];
		r = rho_d[ind];
		fColl_d[ind     ] = omega_d * feqc2D(u,  0, v,  0, r, 4./9.)  + (1.0-omega_d) * f_d[ind     ];
		fColl_d[ind+1*ms] = omega_d * feqc2D(u,  1, v,  0, r, 1./9.)  + (1.0-omega_d) * f_d[ind+1*ms];
		fColl_d[ind+2*ms] = omega_d * feqc2D(u,  0, v,  1, r, 1./9.)  + (1.0-omega_d) * f_d[ind+2*ms];
		fColl_d[ind+3*ms] = omega_d * feqc2D(u, -1, v,  0, r, 1./9.)  + (1.0-omega_d) * f_d[ind+3*ms];
		fColl_d[ind+4*ms] = omega_d * feqc2D(u,  0, v, -1, r, 1./9.)  + (1.0-omega_d) * f_d[ind+4*ms];
		fColl_d[ind+5*ms] = omega_d * feqc2D(u,  1, v,  1, r, 1./36.) + (1.0-omega_d) * f_d[ind+5*ms];
		fColl_d[ind+6*ms] = omega_d * feqc2D(u, -1, v,  1, r, 1./36.) + (1.0-omega_d) * f_d[ind+6*ms];
		fColl_d[ind+7*ms] = omega_d * feqc2D(u, -1, v, -1, r, 1./36.) + (1.0-omega_d) * f_d[ind+7*ms];
		fColl_d[ind+8*ms] = omega_d * feqc2D(u,  1, v, -1, r, 1./36.) + (1.0-omega_d) * f_d[ind+8*ms];
	}
}

__global__ void gpuCollBgkwGC2D(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *r_rho_d, FLOAT_TYPE *b_rho_d, FLOAT_TYPE *u_d,
		FLOAT_TYPE *v_d, FLOAT_TYPE *r_f_d, FLOAT_TYPE *b_f_d, FLOAT_TYPE *r_fColl_d, FLOAT_TYPE *b_fColl_d, int *cg_dir_d){

	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d*length_d;
	int index9;
	FLOAT_TYPE r_r, b_r, r, u, v, chi, omega_temp, color_gradient_x, color_gradient_y;
	FLOAT_TYPE k_r, k_b, k_k, color_gradient_norm, cosin;
	FLOAT_TYPE prod_c_g, r_pert, b_pert;
	FLOAT_TYPE r_CollPert, b_CollPert;
	if (ind < ms)
	{
		u =   u_d[ind];
		v =   v_d[ind];
		r_r = r_rho_d[ind];
		b_r = b_rho_d[ind];
		r = rho_d[ind];

		if (r_omega_d != b_omega_d){
			chi=(r_r - b_r)/r;
			if(chi >= -control_param_d && chi <= control_param_d){
				if (chi > del_d)
					omega_temp=r_omega_d;
				else if (chi <= del_d && chi > 0)
					omega_temp=a1_d + a2_d * chi + a3_d * chi * chi;
				else if (chi <= 0 && chi >= -del_d)
					omega_temp=a1_d + a4_d * chi + a5_d * chi * chi;
				else if (chi < -del_d)
					omega_temp=b_omega_d;
			}
		}
		else
			omega_temp=r_omega_d;

		color_gradient_x = calculateColorGradientX(r_rho_d,b_rho_d, cg_dir_d[ind], ind);
		color_gradient_y = calculateColorGradientY(r_rho_d,b_rho_d, cg_dir_d[ind], ind);


		color_gradient_norm = sqrt(color_gradient_x * color_gradient_x + color_gradient_y * color_gradient_y);
		k_r=r_r/r;
		k_b=b_r/r;
		k_k= beta_d * r_r * b_r/(r * r);


		for (int k=0;k<9;k++){
			if (color_gradient_norm > g_limit_d){
				prod_c_g=cx2D_d[k]*color_gradient_x+cy2D_d[k]*color_gradient_y;
				if (k!=0){
					cosin= prod_c_g / (color_gradient_norm*c_norms_d[k]);
				}
				else
					cosin=0.0;

				// calculate perturbation terms
				r_pert=0.5*r_A_d*color_gradient_norm*(w2D_d[k]*(prod_c_g * prod_c_g) / (color_gradient_norm * color_gradient_norm)-w_pert_d[k]);
				b_pert=0.5*b_A_d*color_gradient_norm*(w2D_d[k]*(prod_c_g * prod_c_g) / (color_gradient_norm * color_gradient_norm)-w_pert_d[k]);
			}
			else{
				// ther perturbation terms are null
				r_pert=0;
				b_pert=0;
			}

			index9 = ind + k * ms;
			// calculate updated distribution function
			if(fluid_d[ind] == 1){
				r_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], r_r, w2D_d[k], r_phi_d[k]) + (1-omega_temp)*r_f_d[index9]+r_pert;
				b_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], b_r, w2D_d[k], b_phi_d[k]) + (1-omega_temp)*b_f_d[index9]+b_pert;
			} else{
				r_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], r_r, w2D_d[k], r_phi_d[k]) + (1-omega_temp)*r_f_d[index9]+r_pert;
				b_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], b_r, w2D_d[k], b_phi_d[k]) + (1-omega_temp)*b_f_d[index9]+b_pert;
			}
			//perform recolor step
			r_fColl_d[index9]=k_r*(r_CollPert + b_CollPert)+k_k*cosin*(r_r*r_phi_d[k]+b_r*b_phi_d[k]);
			b_fColl_d[index9]=k_b*(r_CollPert + b_CollPert)-k_k*cosin*(r_r*r_phi_d[k]+b_r*b_phi_d[k]);
		}
	}

}

__global__ void gpuCollBgkwGC3D(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *r_rho_d, FLOAT_TYPE *b_rho_d, FLOAT_TYPE *u_d,
		FLOAT_TYPE *v_d, FLOAT_TYPE *w_d, FLOAT_TYPE *f_d, FLOAT_TYPE *r_fColl_d, FLOAT_TYPE *b_fColl_d, int *cg_dir_d){

	int ind =  (blockIdx.x + blockIdx.y * gridDim.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
	int ms = depth_d*length_d*height_d;
	int index9;
	FLOAT_TYPE r_r, b_r, r, u, v, w, chi, omega_temp, color_gradient_x, color_gradient_y;
	FLOAT_TYPE k_r, k_b, k_k, color_gradient_norm, cosin;
	FLOAT_TYPE prod_c_g, r_pert, b_pert;
	FLOAT_TYPE r_CollPert, b_CollPert;
	if (ind < ms)
	{
		u =   u_d[ind];
		v =   v_d[ind];
		r_r = r_rho_d[ind];
		b_r = b_rho_d[ind];
		r = rho_d[ind];

		color_gradient_x = calculateColorGradientX(r_rho_d,b_rho_d, cg_dir_d[ind], ind);
		color_gradient_y = calculateColorGradientY(r_rho_d,b_rho_d, cg_dir_d[ind], ind);


		color_gradient_norm = sqrt(color_gradient_x * color_gradient_x + color_gradient_y * color_gradient_y);
		k_r=r_r/r;
		k_b=b_r/r;
		k_k= beta_d * r_r * b_r/(r * r);


		for (int k=0;k<9;k++){
			if (color_gradient_norm > g_limit_d){
				prod_c_g=cx2D_d[k]*color_gradient_x+cy2D_d[k]*color_gradient_y;
				if (k!=0){
					cosin= prod_c_g / (color_gradient_norm*c_norms_d[k]);
				}
				else
					cosin=0.0;

				// calculate perturbation terms
				r_pert=0.5*r_A_d*color_gradient_norm*(w2D_d[k]*(prod_c_g * prod_c_g) / (color_gradient_norm * color_gradient_norm)-w_pert_d[k]);
				b_pert=0.5*b_A_d*color_gradient_norm*(w2D_d[k]*(prod_c_g * prod_c_g) / (color_gradient_norm * color_gradient_norm)-w_pert_d[k]);
			}
			else{
				// ther perturbation terms are null
				r_pert=0;
				b_pert=0;
			}

			index9 = ind + k * ms;
			// calculate updated distribution function
//			if(fluid_d[ind] == 1){
//				r_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], r_r, w2D_d[k], r_phi_d[k]) + (1-omega_temp)*r_f_d[index9]+r_pert;
//				b_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], b_r, w2D_d[k], b_phi_d[k]) + (1-omega_temp)*b_f_d[index9]+b_pert;
//			} else{
//				r_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], r_r, w2D_d[k], r_phi_d[k]) + (1-omega_temp)*r_f_d[index9]+r_pert;
//				b_CollPert = omega_temp*feqc2DCG(u,  cx2D_d[k], v,  cy2D_d[k], b_r, w2D_d[k], b_phi_d[k]) + (1-omega_temp)*b_f_d[index9]+b_pert;
//			}
			//perform recolor step
			r_fColl_d[index9]=k_r*(r_CollPert + b_CollPert)+k_k*cosin*(r_r*r_phi_d[k]+b_r*b_phi_d[k]);
			b_fColl_d[index9]=k_b*(r_CollPert + b_CollPert)-k_k*cosin*(r_r*r_phi_d[k]+b_r*b_phi_d[k]);
		}
	}

}

__global__ void gpuCollBgkw3D(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d,
		FLOAT_TYPE *v_d, FLOAT_TYPE *w_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int ind =  blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;

	int ms = depth_d*length_d*height_d;
	FLOAT_TYPE r, u, v, w;
	if (ind < ms && fluid_d[ind] == 1)
	{
		u = u_d[ind];
		v = v_d[ind];
		w = w_d[ind];
		r = rho_d[ind];

		fColl_d[ ind + 0 *ms ] = omega_d * feqc3D(u, cx3D_d[ 0 ], v, cy3D_d[ 0 ], w, cz3D_d[ 0 ], r, w3D_d[ 0 ]) + (1.0-omega_d) * f_d[ind+ 0 *ms];
		fColl_d[ ind + 1 *ms ] = omega_d * feqc3D(u, cx3D_d[ 1 ], v, cy3D_d[ 1 ], w, cz3D_d[ 1 ], r, w3D_d[ 1 ]) + (1.0-omega_d) * f_d[ind+ 1 *ms];
		fColl_d[ ind + 2 *ms ] = omega_d * feqc3D(u, cx3D_d[ 2 ], v, cy3D_d[ 2 ], w, cz3D_d[ 2 ], r, w3D_d[ 2 ]) + (1.0-omega_d) * f_d[ind+ 2 *ms];
		fColl_d[ ind + 3 *ms ] = omega_d * feqc3D(u, cx3D_d[ 3 ], v, cy3D_d[ 3 ], w, cz3D_d[ 3 ], r, w3D_d[ 3 ]) + (1.0-omega_d) * f_d[ind+ 3 *ms];
		fColl_d[ ind + 4 *ms ] = omega_d * feqc3D(u, cx3D_d[ 4 ], v, cy3D_d[ 4 ], w, cz3D_d[ 4 ], r, w3D_d[ 4 ]) + (1.0-omega_d) * f_d[ind+ 4 *ms];
		fColl_d[ ind + 5 *ms ] = omega_d * feqc3D(u, cx3D_d[ 5 ], v, cy3D_d[ 5 ], w, cz3D_d[ 5 ], r, w3D_d[ 5 ]) + (1.0-omega_d) * f_d[ind+ 5 *ms];
		fColl_d[ ind + 6 *ms ] = omega_d * feqc3D(u, cx3D_d[ 6 ], v, cy3D_d[ 6 ], w, cz3D_d[ 6 ], r, w3D_d[ 6 ]) + (1.0-omega_d) * f_d[ind+ 6 *ms];
		fColl_d[ ind + 7 *ms ] = omega_d * feqc3D(u, cx3D_d[ 7 ], v, cy3D_d[ 7 ], w, cz3D_d[ 7 ], r, w3D_d[ 7 ]) + (1.0-omega_d) * f_d[ind+ 7 *ms];
		fColl_d[ ind + 8 *ms ] = omega_d * feqc3D(u, cx3D_d[ 8 ], v, cy3D_d[ 8 ], w, cz3D_d[ 8 ], r, w3D_d[ 8 ]) + (1.0-omega_d) * f_d[ind+ 8 *ms];
		fColl_d[ ind + 9 *ms ] = omega_d * feqc3D(u, cx3D_d[ 9 ], v, cy3D_d[ 9 ], w, cz3D_d[ 9 ], r, w3D_d[ 9 ]) + (1.0-omega_d) * f_d[ind+ 9 *ms];
		fColl_d[ ind +10 *ms ] = omega_d * feqc3D(u, cx3D_d[10 ], v, cy3D_d[10 ], w, cz3D_d[10 ], r, w3D_d[10 ]) + (1.0-omega_d) * f_d[ind+10 *ms];
		fColl_d[ ind +11 *ms ] = omega_d * feqc3D(u, cx3D_d[11 ], v, cy3D_d[11 ], w, cz3D_d[11 ], r, w3D_d[11 ]) + (1.0-omega_d) * f_d[ind+11 *ms];
		fColl_d[ ind +12 *ms ] = omega_d * feqc3D(u, cx3D_d[12 ], v, cy3D_d[12 ], w, cz3D_d[12 ], r, w3D_d[12 ]) + (1.0-omega_d) * f_d[ind+12 *ms];
		fColl_d[ ind +13 *ms ] = omega_d * feqc3D(u, cx3D_d[13 ], v, cy3D_d[13 ], w, cz3D_d[13 ], r, w3D_d[13 ]) + (1.0-omega_d) * f_d[ind+13 *ms];
		fColl_d[ ind +14 *ms ] = omega_d * feqc3D(u, cx3D_d[14 ], v, cy3D_d[14 ], w, cz3D_d[14 ], r, w3D_d[14 ]) + (1.0-omega_d) * f_d[ind+14 *ms];
		fColl_d[ ind +15 *ms ] = omega_d * feqc3D(u, cx3D_d[15 ], v, cy3D_d[15 ], w, cz3D_d[15 ], r, w3D_d[15 ]) + (1.0-omega_d) * f_d[ind+15 *ms];
		fColl_d[ ind +16 *ms ] = omega_d * feqc3D(u, cx3D_d[16 ], v, cy3D_d[16 ], w, cz3D_d[16 ], r, w3D_d[16 ]) + (1.0-omega_d) * f_d[ind+16 *ms];
		fColl_d[ ind +17 *ms ] = omega_d * feqc3D(u, cx3D_d[17 ], v, cy3D_d[17 ], w, cz3D_d[17 ], r, w3D_d[17 ]) + (1.0-omega_d) * f_d[ind+17 *ms];
		fColl_d[ ind +18 *ms ] = omega_d * feqc3D(u, cx3D_d[18 ], v, cy3D_d[18 ], w, cz3D_d[18 ], r, w3D_d[18 ]) + (1.0-omega_d) * f_d[ind+18 *ms];
	}
}
__global__ void gpuCollTrt(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d,
		FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d*length_d;
	FLOAT_TYPE r, u, v;
	if (ind < ms && fluid_d[ind] == 1)
	{
		u = u_d[ind];
		v = v_d[ind];
		r = rho_d[ind];

		FLOAT_TYPE feq0 = feqc2D(u,  0, v,  0, r, 4./9.);
		FLOAT_TYPE feq1 = feqc2D(u,  1, v,  0, r, 1./9.);
		FLOAT_TYPE feq2 = feqc2D(u,  0, v,  1, r, 1./9.);
		FLOAT_TYPE feq3 = feqc2D(u, -1, v,  0, r, 1./9.);
		FLOAT_TYPE feq4 = feqc2D(u,  0, v, -1, r, 1./9.);
		FLOAT_TYPE feq5 = feqc2D(u,  1, v,  1, r, 1./36.);
		FLOAT_TYPE feq6 = feqc2D(u, -1, v,  1, r, 1./36.);
		FLOAT_TYPE feq7 = feqc2D(u, -1, v, -1, r, 1./36.);
		FLOAT_TYPE feq8 = feqc2D(u,  1, v, -1, r, 1./36.);

		FLOAT_TYPE f0 = f_d[ind];
		FLOAT_TYPE f1 = f_d[ind+1*ms];
		FLOAT_TYPE f2 = f_d[ind+2*ms];
		FLOAT_TYPE f3 = f_d[ind+3*ms];
		FLOAT_TYPE f4 = f_d[ind+4*ms];
		FLOAT_TYPE f5 = f_d[ind+5*ms];
		FLOAT_TYPE f6 = f_d[ind+6*ms];
		FLOAT_TYPE f7 = f_d[ind+7*ms];
		FLOAT_TYPE f8 = f_d[ind+8*ms];

		fColl_d[ind]      = f0 - 0.5 * omega_d * (f0+f0 - feq0-feq0) - 0.5 * omegaA_d * (f0-f0 - feq0+feq0);
		fColl_d[ind+1*ms] = f1 - 0.5 * omega_d * (f1+f3 - feq1-feq3) - 0.5 * omegaA_d * (f1-f3 - feq1+feq3);
		fColl_d[ind+2*ms] = f2 - 0.5 * omega_d * (f2+f4 - feq2-feq4) - 0.5 * omegaA_d * (f2-f4 - feq2+feq4);
		fColl_d[ind+3*ms] = f3 - 0.5 * omega_d * (f3+f1 - feq3-feq1) - 0.5 * omegaA_d * (f3-f1 - feq3+feq1);
		fColl_d[ind+4*ms] = f4 - 0.5 * omega_d * (f4+f2 - feq4-feq2) - 0.5 * omegaA_d * (f4-f2 - feq4+feq2);
		fColl_d[ind+5*ms] = f5 - 0.5 * omega_d * (f5+f7 - feq5-feq7) - 0.5 * omegaA_d * (f5-f7 - feq5+feq7);
		fColl_d[ind+6*ms] = f6 - 0.5 * omega_d * (f6+f8 - feq6-feq8) - 0.5 * omegaA_d * (f6-f8 - feq6+feq8);
		fColl_d[ind+7*ms] = f7 - 0.5 * omega_d * (f7+f5 - feq7-feq5) - 0.5 * omegaA_d * (f7-f5 - feq7+feq5);
		fColl_d[ind+8*ms] = f8 - 0.5 * omega_d * (f8+f6 - feq8-feq6) - 0.5 * omegaA_d * (f8-f6 - feq8+feq6);

		// fColl_d[ind]      = f0 - omega_d * (0.5*(f0+f0) - 0.5*(feq0+feq0)) - omegaA_d * (0.5*(f0-f0) - 0.5*(feq0-feq0));
		// fColl_d[ind+1*ms] = f1 - omega_d * (0.5*(f1+f3) - 0.5*(feq1+feq3)) - omegaA_d * (0.5*(f1-f3) - 0.5*(feq1-feq3));
		// fColl_d[ind+2*ms] = f2 - omega_d * (0.5*(f2+f4) - 0.5*(feq2+feq4)) - omegaA_d * (0.5*(f2-f4) - 0.5*(feq2-feq4));
		// fColl_d[ind+3*ms] = f3 - omega_d * (0.5*(f3+f1) - 0.5*(feq3+feq1)) - omegaA_d * (0.5*(f3-f1) - 0.5*(feq3-feq1));
		// fColl_d[ind+4*ms] = f4 - omega_d * (0.5*(f4+f2) - 0.5*(feq4+feq2)) - omegaA_d * (0.5*(f4-f2) - 0.5*(feq4-feq2));
		// fColl_d[ind+5*ms] = f5 - omega_d * (0.5*(f5+f7) - 0.5*(feq5+feq7)) - omegaA_d * (0.5*(f5-f7) - 0.5*(feq5-feq7));
		// fColl_d[ind+6*ms] = f6 - omega_d * (0.5*(f6+f8) - 0.5*(feq6+feq8)) - omegaA_d * (0.5*(f6-f8) - 0.5*(feq6-feq8));
		// fColl_d[ind+7*ms] = f7 - omega_d * (0.5*(f7+f5) - 0.5*(feq7+feq5)) - omegaA_d * (0.5*(f7-f5) - 0.5*(feq7-feq5));
		// fColl_d[ind+8*ms] = f8 - omega_d * (0.5*(f8+f6) - 0.5*(feq8+feq6)) - omegaA_d * (0.5*(f8-f6) - 0.5*(feq8-feq6));
	}
}

__global__ void gpuCollMrt2D(int* fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d, FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d*length_d;
	FLOAT_TYPE mEq[9], m[9], collision[9], f[9];

	FLOAT_TYPE r,u,v;

	if (ind < ms && fluid_d[ind] == 1)
	{
		r = rho_d[ind];
		u = u_d[ind];
		v = v_d[ind];

		// m^eq = (rho, e, epsilon, jx, qx, jy, qy, pxx, pxy)^T
		mEq[0] = r;
		mEq[1] = r * (-2. + 3. * r * (u*u + v*v));
		mEq[2] = r * (1.  - 3. * r * (u*u + v*v));
		mEq[3] =  r * u;
		mEq[4] = -r * u;
		mEq[5] =  r * v;
		mEq[6] = -r * v;
		mEq[7] = r * (u*u - v*v);
		mEq[8] = r * u * v;

		f[0] = f_d[ind];
		f[1] = f_d[ind+ms];
		f[2] = f_d[ind+2*ms];
		f[3] = f_d[ind+3*ms];
		f[4] = f_d[ind+4*ms];
		f[5] = f_d[ind+5*ms];
		f[6] = f_d[ind+6*ms];
		f[7] = f_d[ind+7*ms];
		f[8] = f_d[ind+8*ms];

		// m = Mf
		for (int i=0; i<9; ++i)
		{
			m[i] = 0;
			for (int j=0; j<9; ++j)
			{
				m[i] += velMomMap2D_d[i*9+j] * f[j];
			}
		}

		// Diff = M^-1 * S * (m - m^eq)
		for (int i=0; i<9; ++i)
		{
			collision[i] = 0;
			for (int j=0; j<9; ++j)
			{
				collision[i] += momCollMtx2D_d[i*9+j] * (m[j] - mEq[j]);
			}
		}

		// fColl = f - M^-1 * S * (m - m^eq) >>> MRT equation
		fColl_d[ind     ] = f[0] - collision[0];
		fColl_d[ind+  ms] = f[1] - collision[1];
		fColl_d[ind+2*ms] = f[2] - collision[2];
		fColl_d[ind+3*ms] = f[3] - collision[3];
		fColl_d[ind+4*ms] = f[4] - collision[4];
		fColl_d[ind+5*ms] = f[5] - collision[5];
		fColl_d[ind+6*ms] = f[6] - collision[6];
		fColl_d[ind+7*ms] = f[7] - collision[7];
		fColl_d[ind+8*ms] = f[8] - collision[8];
	}
}

__global__ void gpuCollMrt3D(int* fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d,
		FLOAT_TYPE *v_d, FLOAT_TYPE *w_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
	int blockId = blockIdx.x
			+ blockIdx.y * gridDim.x;
	int ind =  blockId * (blockDim.x * blockDim.y)
																																										+ (threadIdx.y * blockDim.x)
																																										+ threadIdx.x;

	int ms = depth_d*length_d*height_d;
	FLOAT_TYPE mEq[19], mEq2[19], mEq0[19], m[19], collision[19], f[19];

	FLOAT_TYPE r,u,v,w;
	FLOAT_TYPE jx, jy, jz;
	FLOAT_TYPE wej; // wxx, we

	if (ind < ms && fluid_d[ind] == 1)
	{
		r = rho_d[ind];
		u = u_d[ind];
		v = v_d[ind];
		w = w_d[ind];

		jx = r*u; jy = r*v; jz = r*w;
		wej = -475./63.; //wxx = 0.0;we = 0.0;

		mEq[0 ] = r;
		mEq[1 ] = -11.*r + 19.*(jx*jx + jy*jy + jz*jz);
		mEq[2 ] = wej*(jx*jx + jy*jy + jz*jz);
		mEq[3 ] = jx;
		mEq[4 ] = -(2./3.)*jx;
		mEq[5 ] = jy;
		mEq[6 ] = -(2./3.)*jy;
		mEq[7 ] = jz;
		mEq[8 ] = -(2./3.)*jz;
		mEq[9 ] = (2*jx*jx-(jy*jy+jz*jz));
		mEq[10] = 0.0;
		mEq[11] = jy*jy - jz*jz;
		mEq[12] = 0.0;
		mEq[13] = jx*jy;
		mEq[14] = jy*jz;
		mEq[15] = jx*jz;
		mEq[16] = 0.0;
		mEq[17] = 0.0;
		mEq[18] = 0.0;

		//if(ind==1002)printf("i am MRT\n");

		f[	0	]	=	f_d[	ind+0	*ms	];
		f[	1	]	=	f_d[	ind+1	*ms	];
		f[	2	]	=	f_d[	ind+2	*ms	];
		f[	3	]	=	f_d[	ind+3	*ms	];
		f[	4	]	=	f_d[	ind+4	*ms	];
		f[	5	]	=	f_d[	ind+5	*ms	];
		f[	6	]	=	f_d[	ind+6	*ms	];
		f[	7	]	=	f_d[	ind+7	*ms	];
		f[	8	]	=	f_d[	ind+8	*ms	];
		f[	9	]	=	f_d[	ind+9	*ms	];
		f[	10	]	=	f_d[	ind+10	*ms	];
		f[	11	]	=	f_d[	ind+11	*ms	];
		f[	12	]	=	f_d[	ind+12	*ms	];
		f[	13	]	=	f_d[	ind+13	*ms	];
		f[	14	]	=	f_d[	ind+14	*ms	];
		f[	15	]	=	f_d[	ind+15	*ms	];
		f[	16	]	=	f_d[	ind+16	*ms	];
		f[	17	]	=	f_d[	ind+17	*ms	];
		f[	18	]	=	f_d[	ind+18	*ms	];

		// m = Mf
		for (int i=0; i<19; ++i)
		{
			m[i] = 0.;
			for (int j=0; j<19; ++j)
			{
				m[i] += velMomMap3D_d[i*19+j] * f[j];
			}
		}

		// m = Mf


		for (int i=0; i<19; ++i)
		{
			mEq2[i] = 0.;
			for (int j=0; j<19; ++j)
			{
				mEq2[i] += velMomMap3D_d[i*19+j] * feqc3D(u, cx3D_d[ j ],
						v, cy3D_d[ j ], w, cz3D_d[ j ], r, w3D_d[ j ]);
			}
		}
		//if(ind==0)for (int i=0; i<19; ++i)printf("mEq0(%d)= %.14f\n",i,mEq2[i]-mEq[i]);


		// Diff = M^-1 * S * (m - m^eq)
		for (int i=0; i<19; ++i)
		{
			collision[i] = 0.;
			for (int j=0; j<19; ++j)
			{
				collision[i] += momCollMtx3D_d[i*19+j] * (m[j] - mEq[j]);
			}
		}

		// fColl = f - M^-1 * S * (m - m^eq) >>> MRT equation
		fColl_d[ind  +  0   *ms]  =  f[  0   ]  -  collision[  0   ];
		fColl_d[ind  +  1   *ms]  =  f[  1   ]  -  collision[  1   ];
		fColl_d[ind  +  2   *ms]  =  f[  2   ]  -  collision[  2   ];
		fColl_d[ind  +  3   *ms]  =  f[  3   ]  -  collision[  3   ];
		fColl_d[ind  +  4   *ms]  =  f[  4   ]  -  collision[  4   ];
		fColl_d[ind  +  5   *ms]  =  f[  5   ]  -  collision[  5   ];
		fColl_d[ind  +  6   *ms]  =  f[  6   ]  -  collision[  6   ];
		fColl_d[ind  +  7   *ms]  =  f[  7   ]  -  collision[  7   ];
		fColl_d[ind  +  8   *ms]  =  f[  8   ]  -  collision[  8   ];
		fColl_d[ind  +  9   *ms]  =  f[  9   ]  -  collision[  9   ];
		fColl_d[ind  +  10  *ms]  =  f[  10  ]  -  collision[  10  ];
		fColl_d[ind  +  11  *ms]  =  f[  11  ]  -  collision[  11  ];
		fColl_d[ind  +  12  *ms]  =  f[  12  ]  -  collision[  12  ];
		fColl_d[ind  +  13  *ms]  =  f[  13  ]  -  collision[  13  ];
		fColl_d[ind  +  14  *ms]  =  f[  14  ]  -  collision[  14  ];
		fColl_d[ind  +  15  *ms]  =  f[  15  ]  -  collision[  15  ];
		fColl_d[ind  +  16  *ms]  =  f[  16  ]  -  collision[  16  ];
		fColl_d[ind  +  17  *ms]  =  f[  17  ]  -  collision[  17  ];
		//        fColl_d[ind  +  18  *ms]  =  f[  18  ]  -  collision[  18  ];
		//        if(ind==1123){
		//        	printf("MRT f_coll f0: %.14f f1: %.14f f2: %.14f f3: %.14f f4: %.14f f5: %.14f f6: %.14f f7: %.14f f8: %.14f "
		//					"f9: %.14f f10: %.14f f11: %.14f f12: %.14f f3: %.14f f4: %.14f f15: %.14f f16: %.14f f17: %.14f f18: %.14f\n",
		//					  fColl_d[ind+ 0*ms],fColl_d[ind+ 1*ms],fColl_d[ind+ 2*ms],fColl_d[ind+ 3*ms],fColl_d[ind+ 4*ms],fColl_d[ind+ 5*ms],fColl_d[ind+ 6*ms],fColl_d[ind+ 7*ms], fColl_d[ind+ 8*ms], fColl_d[ind+ 9*ms],
		//					fColl_d[ind+10*ms],fColl_d[ind+11*ms],fColl_d[ind+12*ms],fColl_d[ind+13*ms],fColl_d[ind+14*ms],fColl_d[ind+15*ms], fColl_d[ind+16*ms], fColl_d[ind+17*ms], fColl_d[ind+18*ms]);
		//
		//        }
	}
}

__global__ void gpuCollMrt3D_short(int* fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d,
		FLOAT_TYPE *v_d, FLOAT_TYPE *w_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
	int blockId = blockIdx.x
			+ blockIdx.y * gridDim.x;
	int ind =  blockId * (blockDim.x * blockDim.y)
																																										+ (threadIdx.y * blockDim.x)
																																										+ threadIdx.x;

	int ms = depth_d*length_d*height_d;
	FLOAT_TYPE mEq[19], m[19], collision[19];

	FLOAT_TYPE r,u,v,w;

	if (ind < ms && fluid_d[ind] == 1)
	{
		r = rho_d[ind];
		u = u_d[ind];
		v = v_d[ind];
		w = w_d[ind];


		// m = Mf
		for (int i=0; i<19; ++i)
		{
			m[i] = 0.;
			for (int j=0; j<19; ++j)
			{
				m[i] += velMomMap3D_d[i*19+j] * f_d[ ind + j*ms	];
			}
		}

		// mEq = Mfeq


		for (int i=0; i<19; ++i)
		{
			mEq[i] = 0.;
			for (int j=0; j<19; ++j)
			{
				mEq[i] += velMomMap3D_d[i*19+j] * feqc3D(u, cx3D_d[ j ],
						v, cy3D_d[ j ], w, cz3D_d[ j ], r, w3D_d[ j ]);
			}
		}

		// Diff = M^-1 * S * (m - m^eq)
		for (int i=0; i<19; ++i)
		{
			collision[i] = 0.;
			for (int j=0; j<19; ++j)
			{
				collision[i] += momCollMtx3D_d[i*19+j] * (m[j] - mEq[j]);
				//				if(ind == 0)printf("%f,",momCollMtx3D_d[i*19+j]);
			}
			//			if(ind == 0)printf("\n");
		}

		// fColl = f - M^-1 * S * (m - m^eq) >>> MRT equation
		for (int i=0; i<19; ++i)
		{
			fColl_d[ind  +  i   *ms]  =  f_d[ ind + i*ms	]  -  collision[i];
		}
		if(ind==1123){
			//        	printf("MRT f_coll f0: %.14f f1: %.14f f2: %.14f f3: %.14f f4: %.14f f5: %.14f f6: %.14f f7: %.14f f8: %.14f "
			//					"f9: %.14f f10: %.14f f11: %.14f f12: %.14f f3: %.14f f4: %.14f f15: %.14f f16: %.14f f17: %.14f f18: %.14f\n",
			//					  fColl_d[ind+ 0*ms],fColl_d[ind+ 1*ms],fColl_d[ind+ 2*ms],fColl_d[ind+ 3*ms],fColl_d[ind+ 4*ms],fColl_d[ind+ 5*ms],fColl_d[ind+ 6*ms],fColl_d[ind+ 7*ms], fColl_d[ind+ 8*ms], fColl_d[ind+ 9*ms],
			//					fColl_d[ind+10*ms],fColl_d[ind+11*ms],fColl_d[ind+12*ms],fColl_d[ind+13*ms],fColl_d[ind+14*ms],fColl_d[ind+15*ms], fColl_d[ind+16*ms], fColl_d[ind+17*ms], fColl_d[ind+18*ms]);

		}
	}
}
