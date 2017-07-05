#include <stdio.h>

#include "GpuFunctions.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "GpuConstants.h"

__global__ void gpuStreaming2D(int* fluid_d, int* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d)
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d*length_d;
	FLOAT_TYPE *f, *mf;
	int n = length_d;
	if (ind < ms && fluid_d[ind] == 1)
	{
		f_d[ind] = fColl_d[ind];	//Update fNewStep = fColl
		f = f_d + ms;				// f is f_d memory positions but f starts in f_d 1st level==1st lattice direction
		mf = fColl_d + ms;
//		f[ind]      = (stream_d[ind]      == 1) ? mf[ind-1]        : mf[ind];		// stream_d == 1 means that
//		f[ind+ms]   = (stream_d[ind+ms]   == 1) ? mf[ind+ms-n]     : mf[ind+ms]; 	// the streaming is allowed
//		f[ind+2*ms] = (stream_d[ind+2*ms] == 1) ? mf[ind+2*ms+1]   : mf[ind+2*ms];  // "the regular case"
//		f[ind+3*ms] = (stream_d[ind+3*ms] == 1) ? mf[ind+3*ms+n]   : mf[ind+3*ms];  // stream_d != 1 means
//		f[ind+4*ms] = (stream_d[ind+4*ms] == 1) ? mf[ind+4*ms-n-1] : mf[ind+4*ms]; 	// wall or node outside dom.
//		f[ind+5*ms] = (stream_d[ind+5*ms] == 1) ? mf[ind+5*ms-n+1] : mf[ind+5*ms];
//		f[ind+6*ms] = (stream_d[ind+6*ms] == 1) ? mf[ind+6*ms+n+1] : mf[ind+6*ms];
//		f[ind+7*ms] = (stream_d[ind+7*ms] == 1) ? mf[ind+7*ms+n-1] : mf[ind+7*ms];

		//ASK ANTONIO

		f[ind]      = (stream_d[ind]      == 1) ? mf[ind-1]        : f[ind];		// stream_d == 1 means that
		f[ind+ms]   = (stream_d[ind+ms]   == 1) ? mf[ind+ms-n]     : f[ind+ms]; 	// the streaming is allowed
		f[ind+2*ms] = (stream_d[ind+2*ms] == 1) ? mf[ind+2*ms+1]   : f[ind+2*ms];  // "the regular case"
		f[ind+3*ms] = (stream_d[ind+3*ms] == 1) ? mf[ind+3*ms+n]   : f[ind+3*ms];  // stream_d != 1 means
		f[ind+4*ms] = (stream_d[ind+4*ms] == 1) ? mf[ind+4*ms-n-1] : f[ind+4*ms]; 	// wall or node outside dom.
		f[ind+5*ms] = (stream_d[ind+5*ms] == 1) ? mf[ind+5*ms-n+1] : f[ind+5*ms];
		f[ind+6*ms] = (stream_d[ind+6*ms] == 1) ? mf[ind+6*ms+n+1] : f[ind+6*ms];
		f[ind+7*ms] = (stream_d[ind+7*ms] == 1) ? mf[ind+7*ms+n-1] : f[ind+7*ms];

	}
}

__global__ void gpuStreaming3D(int* fluid_d, bool* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d)
{
	int blockId = blockIdx.x
			+ blockIdx.y * gridDim.x;
	int ind =  blockId * (blockDim.x * blockDim.y)
				+ (threadIdx.y * blockDim.x)
				+ threadIdx.x;

	int ms = depth_d*length_d*height_d;
	FLOAT_TYPE *f, *mf;
	if (ind < ms && fluid_d[ind] == 1)
	{
		f_d[ind] = fColl_d[ind];	//Update fNewStep = fColl
		f = f_d + ms;				// f is f_d memory position but f starts in f_d 1st level==1st lattice direction
		mf = fColl_d + ms;
		f[ind+0  *ms]	=	(stream_d[ind+0	 *ms]	==	1)	?	mf[ind+0  *ms +	c3D_d[1	]]:	mf[ind+0  *ms];
		f[ind+1	 *ms]	=	(stream_d[ind+1	 *ms]	==	1)	?	mf[ind+1  *ms +	c3D_d[2	]]:	mf[ind+1  *ms];
		f[ind+2	 *ms]	=	(stream_d[ind+2	 *ms]	==	1)	?	mf[ind+2  *ms +	c3D_d[3	]]:	mf[ind+2  *ms];
		f[ind+3	 *ms]	=	(stream_d[ind+3	 *ms]	==	1)	?	mf[ind+3  *ms +	c3D_d[4	]]:	mf[ind+3  *ms];
		f[ind+4	 *ms]	=	(stream_d[ind+4	 *ms]	==	1)	?	mf[ind+4  *ms +	c3D_d[5	]]:	mf[ind+4  *ms];
		f[ind+5	 *ms]	=	(stream_d[ind+5	 *ms]	==	1)	?	mf[ind+5  *ms +	c3D_d[6	]]:	mf[ind+5  *ms];
		f[ind+6	 *ms]	=	(stream_d[ind+6	 *ms]	==	1)	?	mf[ind+6  *ms +	c3D_d[7	]]:	mf[ind+6  *ms];
		f[ind+7	 *ms]	=	(stream_d[ind+7	 *ms]	==	1)	?	mf[ind+7  *ms +	c3D_d[8	]]:	mf[ind+7  *ms];
		f[ind+8	 *ms]	=	(stream_d[ind+8	 *ms]	==	1)	?	mf[ind+8  *ms +	c3D_d[9	]]:	mf[ind+8  *ms];
		f[ind+9	 *ms]	=	(stream_d[ind+9	 *ms]	==	1)	?	mf[ind+9  *ms +	c3D_d[10]]:	mf[ind+9  *ms];
		f[ind+10 *ms]	=	(stream_d[ind+10 *ms]	==	1)	?	mf[ind+10 *ms +	c3D_d[11]]:	mf[ind+10 *ms];
		f[ind+11 *ms]	=	(stream_d[ind+11 *ms]	==	1)	?	mf[ind+11 *ms +	c3D_d[12]]:	mf[ind+11 *ms];
		f[ind+12 *ms]	=	(stream_d[ind+12 *ms]	==	1)	?	mf[ind+12 *ms +	c3D_d[13]]:	mf[ind+12 *ms];
		f[ind+13 *ms]	=	(stream_d[ind+13 *ms]	==	1)	?	mf[ind+13 *ms +	c3D_d[14]]:	mf[ind+13 *ms];
		f[ind+14 *ms]	=	(stream_d[ind+14 *ms]	==	1)	?	mf[ind+14 *ms +	c3D_d[15]]:	mf[ind+14 *ms];
		f[ind+15 *ms]	=	(stream_d[ind+15 *ms]	==	1)	?	mf[ind+15 *ms +	c3D_d[16]]:	mf[ind+15 *ms];
		f[ind+16 *ms]	=	(stream_d[ind+16 *ms]	==	1)	?	mf[ind+16 *ms +	c3D_d[17]]:	mf[ind+16 *ms];
		f[ind+17 *ms]	=	(stream_d[ind+17 *ms]	==	1)	?	mf[ind+17 *ms +	c3D_d[18]]:	mf[ind+17 *ms];
	}
}
