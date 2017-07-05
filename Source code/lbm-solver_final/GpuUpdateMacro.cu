#include "GpuFunctions.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "GpuConstants.h"

__global__ void gpuUpdateMacro2D(int *fluid_d, FLOAT_TYPE* rho_d,
		FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, int *bcMask_d, FLOAT_TYPE* drag_d,
		FLOAT_TYPE* lift_d,
		FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d) {
	int ind = threadIdx.x + blockIdx.x * blockDim.x;

	int ms = depth_d * length_d;

	FLOAT_TYPE r, u, v;

	if (ind < ms) {
		if (fluid_d[ind] == 1) {
			r = u = v = 0.0;
			r = f_d[ind] + f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
			                                                       + f_d[ind + 4 * ms] + f_d[ind + 5 * ms] + f_d[ind + 6 * ms]
			                                                                                                     + f_d[ind + 7 * ms] + f_d[ind + 8 * ms];
			u = f_d[ind + ms] - f_d[ind + 3 * ms] + f_d[ind + 5 * ms]
			                                            - f_d[ind + 6 * ms] - f_d[ind + 7 * ms] + f_d[ind + 8 * ms];
			v = f_d[ind + 2 * ms] - f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
			                                                + f_d[ind + 6 * ms] - f_d[ind + 7 * ms] - f_d[ind + 8 * ms];

			rho_d[ind] = r;
			u_d[ind] = u / r;
			///@todo code: probably should handle outlet on other sides
			v_d[ind] = ((bcMask_d[ind] & BC_OUTL_E) == BC_OUTL_E) ? 0.0 : v / r;

			//   DRAG/LIFT FORCE
			if (dlBoundaryId_d
					!= 0&& (bcMask_d[ind] & BND_ID_ALL) == BOUND_ID(dlBoundaryId_d)) {
				// printf("draglift: %d\n",ind);
				drag_d[ind] = 0.33333333 * r * (20 - coordX_d[ind]) * 0.2;
				lift_d[ind] = 0.33333333 * r * (20 - coordY_d[ind]) * 0.2;
			}
		}
	}
}

__global__ void gpuUpdateMacro2DCG(int *fluid_d, FLOAT_TYPE* rho_d,
		FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* r_f_d, FLOAT_TYPE* b_f_d, FLOAT_TYPE* r_rho_d,
		FLOAT_TYPE* b_rho_d, FLOAT_TYPE *p_in_d, FLOAT_TYPE *p_out_d,
		int *num_in_d, int *num_out_d) {
	int ind = threadIdx.x + blockIdx.x * blockDim.x;

	int ms = depth_d * length_d;

	FLOAT_TYPE r_r, b_r, u, v, r, chi;

	if (ind < ms) {
		//necessary because of sum
		p_in_d[ind] = 0;
		p_out_d[ind] = 0;
		num_in_d[ind] = 0;
		num_out_d[ind] = 0;

		if (fluid_d[ind] == 1) {
			r_r = b_r = u = v = 0.0;

			r_r = r_f_d[ind] +
					r_f_d[ind + ms] +
					r_f_d[ind + 2 * ms] +
					r_f_d[ind + 3 * ms] +
					r_f_d[ind + 4 * ms] +
					r_f_d[ind + 5 * ms] +
					r_f_d[ind + 6 * ms]	+
					r_f_d[ind + 7 * ms] +
					r_f_d[ind + 8 * ms];
			b_r = b_f_d[ind] +
					b_f_d[ind + ms] +
					b_f_d[ind + 2 * ms] +
					b_f_d[ind + 3 * ms] +
					b_f_d[ind + 4 * ms] +
					b_f_d[ind + 5 * ms] +
					b_f_d[ind + 6 * ms] +
					b_f_d[ind + 7 * ms] +
					b_f_d[ind + 8 * ms];

			r_rho_d[ind] = r_r;
			b_rho_d[ind] = b_r;
			r = r_r + b_r;
			rho_d[ind] = r;

			u = (r_f_d[ind + ms] + b_f_d[ind + ms]) -
					(r_f_d[ind + 3 * ms] + b_f_d[ind + 3 * ms]) +
					(r_f_d[ind + 5 * ms] + b_f_d[ind + 5 * ms]) -
					(r_f_d[ind + 6 * ms] + b_f_d[ind + 6 * ms]) -
					(r_f_d[ind + 7 * ms] + b_f_d[ind + 7 * ms]) +
					(r_f_d[ind + 8 * ms] + b_f_d[ind + 8 * ms]);

			v = (r_f_d[ind + 2 * ms] + b_f_d[ind + 2 * ms]) -
					(r_f_d[ind + 4 * ms] + b_f_d[ind + 4 * ms]) +
					(r_f_d[ind + 5 * ms] + b_f_d[ind + 5 * ms]) +
					(r_f_d[ind + 6 * ms] + b_f_d[ind + 6 * ms]) -
					(r_f_d[ind + 7 * ms] + b_f_d[ind + 7 * ms]) -
					(r_f_d[ind + 8 * ms] + b_f_d[ind + 8 * ms]);

			u_d[ind] = u / r;

			v_d[ind] = v / r;

			// p_in and p_out for the surface tension
			chi=(r_r-b_r)/r;

			if (chi >= control_param_d){
				num_in_d[ind] = 1;
				p_in_d[ind] = r_r;
			}
			else if (chi <= -control_param_d){
				num_out_d[ind] = 1;
				p_out_d[ind] = b_r;
			}
		}
	}
}

__global__ void gpuUpdateMacro3D(int *fluid_d, FLOAT_TYPE* rho_d,
		FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* w_d, int* bcBoundId_d,
		FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* coordZ_d,
		FLOAT_TYPE* f_d, FLOAT_TYPE g, unsigned long long *bcMask_d,int updateInltOutl) //	 FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d, FLOAT_TYPE* latF_d,
{
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int ind = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
															+ threadIdx.x;
	int ms = depth_d * length_d * height_d;

	FLOAT_TYPE r, rU, rV, rW;

	if (ind < ms) {

		//    printf("bcMask[ind] |= BC3D_MASK((unsigned long long)bcType[bci], dir); %#016lX\n", (bcMask_d[bci] & BC3D_FLUID) );
		if (fluid_d[ind] == 1
				&& (!(((bcMask_d[ind] & BC3D_OUTL_B) ==BC3D_INLT_B))||updateInltOutl)) {
			//((bcMask_d[ind] & BC3D_OUTL_B) == BC3D_OUTL_B)||
			r = rU = rV = rW = 0.0;
			r = f_d[ind] + f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
			                                                       + f_d[ind + 4 * ms] + f_d[ind + 5 * ms] + f_d[ind + 6 * ms]
			                                                                                                     + f_d[ind + 7 * ms] + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
			                                                                                                                                                   + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
			                                                                                                                                                                              + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
			                                                                                                                                                                                                         + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
			                                                                                                                                                                                                                                    + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
			                                                                                                                                                                                                                                                               + f_d[ind + 18 * ms];
			//                        if(ind==1) printf("i am here BOTTOM WEST_1 r= %.14f f0: %.14f f1: %.14f f2: %.14f f3: %.14f f4: %.14f f5: %.14f f6: %.14f f7: %.14f f8: %.14f "
			//                        		"f9: %.14f f10: %.14f f11: %.14f f12: %.14f f3: %.14f f4: %.14f f15: %.14f f16: %.14f f17: %.14f f18: %.14f\n",
			//                        		r, f_d[ind+ 0*ms],f_d[ind+ 1*ms],f_d[ind+ 2*ms],f_d[ind+ 3*ms],f_d[ind+ 4*ms],f_d[ind+ 5*ms],f_d[ind+ 6*ms],f_d[ind+ 7*ms], f_d[ind+ 8*ms], f_d[ind+ 9*ms],
			//                        		f_d[ind+10*ms],f_d[ind+11*ms],f_d[ind+12*ms],f_d[ind+13*ms],f_d[ind+14*ms],f_d[ind+15*ms], f_d[ind+16*ms], f_d[ind+17*ms], f_d[ind+18*ms]);

			rU = f_d[ind + ms] - f_d[ind + 2 * ms] + f_d[ind + 7 * ms]
			                                             - f_d[ind + 8 * ms] + f_d[ind + 9 * ms] - f_d[ind + 10 * ms]
			                                                                                           + f_d[ind + 11 * ms] - f_d[ind + 12 * ms]
			                                                                                                                      + f_d[ind + 13 * ms] - f_d[ind + 14 * ms];
			//            if(ind==1) printf("i am here TOP lid u= %.14f f1: %.14f f2: %.14f f7: %.14f f8: %.14f "
			//                        		"f9: %.14f f10: %.14f f11: %.14f f12: %.14f f3: %.14f f4: %.14f  \n",
			//                        		rU/r, f_d[ind+ 1*ms],f_d[ind+ 2*ms],f_d[ind+ 7*ms], f_d[ind+ 8*ms], f_d[ind+ 9*ms],
			//                        		f_d[ind+10*ms],f_d[ind+11*ms],f_d[ind+12*ms],f_d[ind+13*ms],f_d[ind+14*ms]);
			rV = f_d[ind + 3 * ms] - f_d[ind + 4 * ms] + f_d[ind + 7 * ms]
			                                                 + f_d[ind + 8 * ms] - f_d[ind + 9 * ms] - f_d[ind + 10 * ms]
			                                                                                               + f_d[ind + 15 * ms] - f_d[ind + 16 * ms]
			                                                                                                                          + f_d[ind + 17 * ms] - f_d[ind + 18 * ms];
			//
			//            			if(ind==30) printf("i am here BOTTOM WEST_1 v= %.14f f3: %.14f f4: %.14f f7: %.14f f8: %.14f f9: %.14f f10: %.14f f15: %.14f f16: %.14f f17: %.14f f18: %.14f\n",rV/r, f_d[ind+ 3*ms],f_d[ind+ 4*ms],f_d[ind+ 7*ms], f_d[ind+ 8*ms], f_d[ind+ 9*ms], f_d[ind+10*ms],
			//            	            	 f_d[ind+15*ms], f_d[ind+16*ms], f_d[ind+17*ms], f_d[ind+18*ms]);
			rW = f_d[ind + 5 * ms] - f_d[ind + 6 * ms] + f_d[ind + 11 * ms]
			                                                 + f_d[ind + 12 * ms] - f_d[ind + 13 * ms]
			                                                                            - f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
			                                                                                                       + f_d[ind + 16 * ms] - f_d[ind + 17 * ms]
			                                                                                                                                  - f_d[ind + 18 * ms];
			rho_d[ind] = r;
			u_d[ind] = rU / r + g / (omega_d);
			;
			v_d[ind] = rV / r;
			w_d[ind] = rW / r;
			///@todo code: probably should handle outlet on other sides
			//			v_d[ind] =
			//					((bcMask_d[ind] & BC3D_OUTL_1) == BC3D_OUTL_1) ?
			//							0.0 : rV / r;
			//			w_d[ind] =
			//					((bcMask_d[ind] & BC3D_OUTL_1) == BC3D_OUTL_1) ?
			//							0.0 : rW / r;
			//            if(ind==1 )printf(" macro %d %f   %f   %f\n",ind, u_d[ind],v_d[ind],w_d[ind]);

			//            if(u_d[ind]>0.0 | u_d[ind]<0.0) {printf("ind: %d, x: %lf y %lf z %lf\n", ind, coordX_d[ind],coordY_d[ind],coordZ_d[ind]); return ;}
			//			if(u_d[ind]>0.0 | u_d[ind]<0.0) {printf("ind: %d, x: %lf y %lf z %.14f\n", ind, coordX_d[ind],coordY_d[ind],coordZ_d[ind]); return ;}
			//			if(u_d[ind]>0.0 | u_d[ind]<0.0) {printf("ind: %d, x: %lf y %lf z %lf\n", ind, coordX_d[ind],coordY_d[ind],coordZ_d[ind]); return ;}

			//   DRAG/LIFT/LATERAL FORCES TODO: find reference and check

			//			if (dlBoundaryId_d != 0 && bcBoundId_d[ind] == dlBoundaryId_d) {
			// printf("draglift: %d\n",ind);
			//                drag_d[ind] = 0.33333333*r*(20-coordX_d[ind])*0.2;
			//                lift_d[ind] = 0.33333333*r*(20-coordZ_d[ind])*0.2;
			//                latF_d[ind] = 0.33333333*r*(20-coordY_d[ind])*0.2;
			//			}
		}
	}

}
//            if(ind==1) printf("i am here BOTTOM WEST_1 r= %.14f f0: %.14f f1: %.14f f2: %.14f f3: %.14f f4: %.14f f5: %.14f f6: %.14f f7: %.14f f8: %.14f "
//            		"f9: %.14f f10: %.14f f11: %.14f f12: %.14f f3: %.14f f4: %.14f f15: %.14f f16: %.14f f17: %.14f f18: %.14f\n",
//            		r, f_d[ind+ 0*ms],f_d[ind+ 1*ms],f_d[ind+ 2*ms],f_d[ind+ 3*ms],f_d[ind+ 4*ms],f_d[ind+ 5*ms],f_d[ind+ 6*ms],f_d[ind+ 7*ms], f_d[ind+ 8*ms], f_d[ind+ 9*ms],
//            		f_d[ind+10*ms],f_d[ind+11*ms],f_d[ind+12*ms],f_d[ind+13*ms],f_d[ind+14*ms],f_d[ind+15*ms], f_d[ind+16*ms], f_d[ind+17*ms], f_d[ind+18*ms]);
//			if(ind==30) printf("i am here BOTTOM WEST_1 v= %.14f f3: %.14f f4: %.14f f7: %.14f f8: %.14f f9: %.14f f10: %.14f f15: %.14f f16: %.14f f17: %.14f f18: %.14f\n",rV/r, f_d[ind+ 3*ms],f_d[ind+ 4*ms],f_d[ind+ 7*ms], f_d[ind+ 8*ms], f_d[ind+ 9*ms], f_d[ind+10*ms],
//	            	 f_d[ind+15*ms], f_d[ind+16*ms], f_d[ind+17*ms], f_d[ind+18*ms]); 26739
//           if(ind==3) printf("i am here BOTTOM South_1 rW= %.14f f5: %.14f f6: %.14f f11: %.14f f12: %.14f f13: %.14f f14: %.14f f15: %.14f "
//            					"f16: %.14f f17: %.14f f18: %.14f\n",rW, f_d[ind+ 5*ms], f_d[ind+ 6*ms], f_d[ind+11*ms], f_d[ind+12*ms], f_d[ind+13*ms], f_d[ind+14*ms],
//            	            	 f_d[ind+15*ms], f_d[ind+16*ms], f_d[ind+17*ms], f_d[ind+18*ms] );
