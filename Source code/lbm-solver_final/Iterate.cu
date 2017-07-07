#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <errno.h>
#include "GpuFunctions.h"       // GPU kernels
#include "ShellFunctions.h"     // For convenience
#include "FilesReading.h"       // For reading files
#include "FilesWriting.h"       // For writing files e.g. tecplot
#include "CellFunctions.h"      // For cell modifications
#include "ComputeResiduals.h"   // residuals
#include "LogWriter.h"
#include "Iterate.h"
#include "ArrayUtils.h"
#include "Check.h"
#include "Multiphase.h"
#include "GpuSum.h"

#define CUDA 1

int Iterate2D(InputFilenames *inFn, Arguments *args) {
	// Time measurement: declaration, begin
	clock_t tStart = clock();

	FILE* logFile;               // file for log
	char autosaveFilename[768];  // autosave filename
	char outputFilename[768];    // initial data will be written to this file
	char finalFilename[768];     // final data will be written to this file
	char logFilename[768];       // path of the .log file
	char residualsFilename[768]; // path of the residuals file
	char timeFilename[768];      // path of time measurement file

	logFilename[0] = '\0';
	residualsFilename[0] = '\0';
	timeFilename[0] = '\0';

	if (strlen(inFn->result)) {
		strcat(logFilename, inFn->result);
		strcat(residualsFilename, inFn->result);
		strcat(timeFilename, inFn->result);
	}
	strcat(logFilename, "lbmsolver.log");
	strcat(residualsFilename, "residuals.dat");
	strcat(timeFilename, "runtimes.dat");

	int autosaveIt = 1; // autosave i variable, will be incremented after every autosave
	int numNodes, numConns; // This will store the number of lines of the read files
	FLOAT_TYPE delta;          // grid spacing
	int n, m;                   // number of nodes in the x and y directions
	FLOAT_TYPE maxInletCoordY; // maximum inlet coordinate in y
	FLOAT_TYPE minInletCoordY; // minimum inlet coordinate in y
	int numInletNodes;         // number of inlet nodes

	int *nodeIdX, *nodeIdY, *nodeType, *bcNodeIdX, *bcNodeIdY, *latticeId,
	*bcType, *bcBoundId,*tempi;
	FLOAT_TYPE *nodeX, *nodeY, *bcX, *bcY,*temp;

	FLOAT_TYPE taskTime[9];
	int i;
	for (i = 0; i < 9; ++i) {
		taskTime[i] = 0.0;
	}

	clock_t tInstant1, tInstant2; // Time measurement points, universal
	clock_t tIterStart, tIterEnd; // Time measurement points: main loop

	// cuda time measurement variables
	cudaEvent_t start, stop;
	float cudatime;
	CHECK(cudaEventCreate(&start));
	CHECK(cudaEventCreate(&stop));
	numNodes = readNodeFile(inFn->node, &nodeIdX, &nodeIdY, &tempi, &nodeX,
			&nodeY, &temp, &nodeType,args->TypeOfProblem);

	if (numNodes == 0) {
		printf("NODES NOT FOUND in file\n");
		return 2;
	}
	int *fluid_d = createGpuArrayInt(numNodes, ARRAY_COPY, 0, nodeType);
	FLOAT_TYPE *coordX_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeX);
	FLOAT_TYPE *coordY_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeY);
	numConns = readConnFile(inFn->bc, &bcNodeIdX, &bcNodeIdY, &tempi,
			&latticeId, &bcType, &bcX, &bcY, &temp, &bcBoundId,args->TypeOfProblem);

	if (numConns == 0) {
		printf("NEIGHBOURING NOT FOUND in file\n");
		return 2;
	}

	int *bcNodeIdX_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdX);
	int *bcNodeIdY_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdX); //WHATCH OUT IdX???
	int *latticeId_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, latticeId);
	int *bcType_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcType);
	int *bcBoundId_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcBoundId);
	FLOAT_TYPE *bcX_d = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcX);
	FLOAT_TYPE *bcY_d = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcY);
	m = getLastValue(nodeIdY, numNodes);
	n = getLastValue(nodeIdX, numNodes);


	delta = getGridSpacing(nodeIdX, nodeIdY, nodeX, numNodes);
	numInletNodes = getNumInletNodes(bcType, latticeId, numConns,
			args->TypeOfProblem);
	maxInletCoordY = getMaxInletCoordY(bcType, latticeId, bcY, delta, numConns,args->TypeOfProblem);
	minInletCoordY = getMinInletCoordY(bcType, latticeId, bcY, delta, numConns,args->TypeOfProblem);
	FLOAT_TYPE *nodeZ = createHostArrayFlt(m * n, ARRAY_ZERO);
	writeInitLog(logFilename, args, delta, m, n, 1, numInletNodes,
			maxInletCoordY, minInletCoordY, 0.0, 0.0);
	logFile = fopen(logFilename, "a");

	// In case of no autosave
	sprintf(autosaveFilename, "NOWHERE!");
	initConstants2D(args, maxInletCoordY, minInletCoordY, delta, m, n);

	dim3 tpb(THREADS); 					 // THREADS/block
	dim3 bpg1((int) (m * n / THREADS) + 1);     // blocks/grid  MxN
	dim3 bpg8((int) (8 * m * n / THREADS) + 1);     // blocks/grid 8MxN
	dim3 bpg9((int) (9 * m * n / THREADS) + 1);     // blocks/grid 9MxN
	dim3 bpgBC((int) (numConns / THREADS) + 1); // blocks/grid N_BC
	// residuals
	FLOAT_TYPE *norm = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *dragSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *liftSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);

	fprintf(logFile, "\n:::: Initializing ::::\n");
	printf("\n:::: Initializing ::::\n");
	CHECK(cudaEventRecord(start, 0));
	FLOAT_TYPE *u = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *v = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *w = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *rho = createHostArrayFlt(m * n, ARRAY_ZERO);

	//Multiphase
	FLOAT_TYPE *r_rho = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *b_rho = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *st_error = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *color_gradient = createHostArrayFlt(m * n * 2, ARRAY_ZERO);
	FLOAT_TYPE *r_f = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_f = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *r_fColl = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_fColl = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE r_omega = 1.0/(3.0 * args->r_viscosity+0.5);
	FLOAT_TYPE	b_omega = 1.0/(3.0 * args->b_viscosity+0.5);
	FLOAT_TYPE st_predicted = (2.0/9.0)*(1.0+1.0/args->gamma)/(0.5*(r_omega+b_omega))*0.5*args->r_density*(args->r_A+args->b_A);
	int *cg_directions = createHostArrayInt(n*m, ARRAY_ZERO);

#if !CUDA
	FLOAT_TYPE r_phi[9];
	FLOAT_TYPE b_phi[9];
	FLOAT_TYPE w_pert[9];
	int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	FLOAT_TYPE weight[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
	if(args->multiPhase){



		int i;
		w_pert[0] = -4.0/ 27.0;
		for(i = 1; i < 5; i++)
			w_pert[i] = 2.0 / 27.0;
		for(i = 5; i < 9; i++)
			w_pert[i] = 5.0 / 108.0;

		r_phi[0] = args->r_alpha;
		for(i = 1; i < 5; i++)
			r_phi[i] = (1.0 - args->r_alpha) / 5.0;
		for(i = 5; i < 9; i++)
			r_phi[i] = (1.0 - args->r_alpha) / 20.0;

		b_phi[0] = args->b_alpha;
		for(i = 1; i < 5; i++)
			b_phi[i] = (1.0 - args->b_alpha) / 5.0;
		for(i = 5; i < 9; i++)
			b_phi[i] = (1.0 - args->b_alpha) / 20.0;


		createBubble(nodeX, nodeY,n,m,args->bubble_radius, r_f, b_f,r_rho,b_rho, args->r_density, args->b_density, r_phi, b_phi, rho);
	}
#endif
	FLOAT_TYPE *rho_d = createGpuArrayFlt(m * n, ARRAY_FILL, args->rho);
	FLOAT_TYPE *r_rho_d = createGpuArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *b_rho_d = createGpuArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *r_f_d = createGpuArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_f_d = createGpuArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *r_fColl_d = createGpuArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_fColl_d = createGpuArrayFlt(m * n * 9, ARRAY_ZERO);
	int *cg_dir_d = createGpuArrayInt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *p_in_d = createGpuArrayFlt(n*m, ARRAY_ZERO);
	FLOAT_TYPE *p_out_d = createGpuArrayFlt(n*m, ARRAY_ZERO);
	FLOAT_TYPE p_in_mean;
	FLOAT_TYPE p_out_mean;
	FLOAT_TYPE ms = n * m;
	int *num_in_d = createGpuArrayInt(n*m, ARRAY_ZERO);
	int *num_out_d = createGpuArrayInt(n*m, ARRAY_ZERO);

	FLOAT_TYPE *u0_d, *v0_d;

	if (args->inletProfile == NO_INLET) {
		u0_d = createGpuArrayFlt(m * n, ARRAY_FILL, args->u);
		v0_d = createGpuArrayFlt(m * n, ARRAY_FILL, args->v);
	} else {
		u0_d = createGpuArrayFlt(m * n, ARRAY_ZERO);
		v0_d = createGpuArrayFlt(m * n, ARRAY_ZERO);
	}
	if (args->inletProfile == INLET) {
		gpuInitInletProfile2D<<<bpg1, tpb>>>(u0_d, v0_d, coordY_d, m * n);
	}
	FLOAT_TYPE *drag_d = createGpuArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *lift_d = createGpuArrayFlt(m * n, ARRAY_ZERO);

	FLOAT_TYPE *f_d = createGpuArrayFlt(9 * m * n, ARRAY_ZERO);
	FLOAT_TYPE *fColl_d = createGpuArrayFlt(9 * m * n, ARRAY_ZERO);

	FLOAT_TYPE *temp9a_d = createGpuArrayFlt(9 * m * n, ARRAY_ZERO);
	FLOAT_TYPE *temp9b_d = createGpuArrayFlt(9 * m * n, ARRAY_ZERO);
	FLOAT_TYPE *tempA_d = createGpuArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *tempB_d = createGpuArrayFlt(m * n, ARRAY_ZERO);

#if CUDA
	if(args->multiPhase){
		initColorGradient(cg_directions, n, m);
		CHECK(cudaMemcpy(cg_dir_d, cg_directions, SIZEINT(m*n), cudaMemcpyHostToDevice));
		initCGBubble<<<bpg1,tpb>>>(coordX_d,coordY_d,r_rho_d, b_rho_d, rho_d, r_f_d, b_f_d);
	}
#endif
	int *mask = createHostArrayInt(m * n, ARRAY_ZERO);
	int *bcMask = createHostArrayInt(m * n, ARRAY_ZERO);
	int *bcIdx = createHostArrayInt(m * n, ARRAY_ZERO);

	FLOAT_TYPE *u_d = createGpuArrayFlt(m * n, ARRAY_CPYD, 0, u0_d);
	FLOAT_TYPE *v_d = createGpuArrayFlt(m * n, ARRAY_CPYD, 0, v0_d);
	int *stream = createHostArrayInt(8 * m * n, ARRAY_FILL, 1);
	FLOAT_TYPE *q = createHostArrayFlt(8 * m * n, ARRAY_FILL, 0.5);

	int bcCount = initBoundaryConditions2D(bcNodeIdX, bcNodeIdY, q, bcBoundId,
			nodeType, bcX, bcY, nodeX, nodeY, latticeId, stream, bcType, bcMask,
			bcIdx, mask, delta, m, n, numConns);

	int *bcIdxCollapsed_d = createGpuArrayInt(bcCount, ARRAY_ZERO);
	int *bcMaskCollapsed_d = createGpuArrayInt(bcCount, ARRAY_ZERO);
	FLOAT_TYPE *qCollapsed_d = createGpuArrayFlt(8 * bcCount, ARRAY_ZERO);

	dim3 bpgB((int) (bcCount / THREADS) + 1); // blocks/grid
	int *bcMask_d = createGpuArrayInt(m * n, ARRAY_COPY, 0, bcMask);
	int *bcIdx_d = createGpuArrayInt(m * n, ARRAY_COPY, 0, bcIdx);

	collapseBc2D(bcIdx, bcIdxCollapsed_d, bcMask, bcMaskCollapsed_d, q,
			qCollapsed_d, mask, m, n, bcCount);

	int *stream_d = createGpuArrayInt(8 * m * n, ARRAY_COPY, 0, stream);


//	CHECK(cudaMemcpy(u, u_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
//	CHECK(cudaMemcpy(v, v_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
//	CHECK(cudaMemcpy(rho, rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));

	CHECK(cudaEventRecord(stop, 0));
	CHECK(cudaEventSynchronize(stop));
	CHECK(cudaEventElapsedTime(&cudatime, start, stop));
	taskTime[T_INIT] += cudatime / 1000;
	fclose(logFile);
	writeNodeNumbers(logFilename, numNodes, numConns, bcCount);
	logFile = fopen(logFilename, "a");


	void *hostArrays[] = { nodeIdX, nodeIdY, nodeX, nodeY, nodeType, bcNodeIdX,
			bcNodeIdY, latticeId, bcType, bcX, bcY, bcBoundId, u, v, rho, mask,
			bcMask, bcIdx, stream, q, norm, dragSum, liftSum, r_rho, b_rho,
			color_gradient,r_f,b_f,r_fColl, b_fColl, st_error, cg_directions};
	void *gpuArrays[] = { coordX_d, coordY_d, fluid_d, bcNodeIdX_d, bcNodeIdY_d,
			latticeId_d, bcType_d, bcX_d, bcY_d, bcBoundId_d, u_d, v_d, rho_d,
			u0_d, v0_d, drag_d, lift_d, f_d, fColl_d, temp9a_d, temp9b_d,
			tempA_d, tempB_d, bcMask_d, bcMaskCollapsed_d, bcIdx_d,
			bcIdxCollapsed_d, stream_d,  qCollapsed_d, r_f_d, r_fColl_d, b_f_d,
			b_fColl_d, cg_dir_d, r_rho_d, b_rho_d, p_in_d, p_out_d, num_in_d, num_out_d};
	fprintf(logFile, "\n:::: Initialization done! ::::\n");

	printf("Initialization took %f seconds\n", taskTime[T_INIT]);

	// Write Initialized data
	switch (args->outputFormat) {
	case CSV:
		sprintf(finalFilename, "%sFinalData.csv", inFn->result);
		break;
	case TECPLOT:
		sprintf(finalFilename, "%sFinalData.dat", inFn->result);
		break;
	case PARAVIEW:
		sprintf(finalFilename, "%sFinalData.vti", inFn->result);
		break;
	}

	tInstant1 = clock(); // Start measuring time
	if(args->multiPhase){
		WriteResultsMultiPhase(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,r_rho,b_rho, nodeType,
				n, m, 1, args->outputFormat);
	}
	else
		WriteResults3D(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho, nodeType,
				n, m, 1, args->outputFormat);
	tInstant2 = clock();
	taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1) / CLOCKS_PER_SEC;
	printf("\nInitialized data was written to %s\n", outputFilename);

	////////////////// ITERATION ///////////////////////

	fprintf(logFile, "\n:::: Start Iterations ::::\n");
	printf("\n:::: Start Iterations ::::\n");

	printf("%d is the number of iterations \n", args->iterations);

	tIterStart = clock(); // Start measuring time of main loop
	int iter = 0;
	while (iter < args->iterations) {


		//		if(args->multiPhase){
		//			CHECK(cudaMemcpy(r_rho_d, r_rho, SIZEFLT(m * n), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(b_rho_d, b_rho, SIZEFLT(m * n), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(r_f_d, r_f, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(b_f_d, b_f, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(r_fColl_d, r_fColl, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(b_fColl_d, b_fColl, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(u_d, u, SIZEFLT(m * n), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(v_d, v, SIZEFLT(m * n), cudaMemcpyHostToDevice));
		//			CHECK(cudaMemcpy(rho_d, rho, SIZEFLT(m * n), cudaMemcpyHostToDevice));
		//		}
		CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0)); // Start measuring time
		switch (args->collisionModel) {

		case BGKW:
			if(args->multiPhase){
				//Collision

#if !CUDA
				mp2DColl(n, m, rho, u, v, r_f, b_f, r_rho, b_rho, r_phi, b_phi, w_pert, color_gradient,
						r_omega, b_omega, args->control_param, args->del, args->beta,
						args->g_limit, args->r_A, args->b_A, r_fColl, b_fColl, weight, cx, cy);
#else
				gpuCollBgkwGC2D<<<bpg1, tpb>>>(fluid_d, rho_d, r_rho_d, b_rho_d, u_d, v_d, r_f_d, b_f_d, r_fColl_d, b_fColl_d, cg_dir_d);
#endif
			}else{
				gpuCollBgkw2D<<<bpg1, tpb>>>(fluid_d, rho_d, u_d, v_d, f_d,
						fColl_d);
			}
			break;
		case TRT:
			gpuCollTrt<<<bpg1, tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
			break;

		case MRT:
			gpuCollMrt2D<<<bpg1, tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
			break;
		}

		CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_COLL] += cudatime;

		////////////// STREAMING ///////////////
		CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0));
		if(args->multiPhase){
#if !CUDA
			streamMP(n, m, r_f, b_f, r_fColl, b_fColl);
			//			gpuStreaming2D<<<bpg1, tpb>>>(fluid_d, stream_d, r_f_d, r_fColl_d);
			//			gpuStreaming2D<<<bpg1, tpb>>>(fluid_d, stream_d, b_f_d, b_fColl_d);
#else
			gpuStreaming2DCG<<<bpg1, tpb>>>(fluid_d, stream_d, r_f_d, r_fColl_d, b_f_d, b_fColl_d);
#endif
		}
		else{
			gpuStreaming2D<<<bpg1, tpb>>>(fluid_d, stream_d, f_d, fColl_d);
		}

		CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_STRM] += cudatime;


		// make the host block until the device is finished with foo
		CHECK(cudaThreadSynchronize());

		// check for error
		cudaError_t error = cudaGetLastError();
		if (error != cudaSuccess) {
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			exit(-1);
		}
		////////////// BOUNDARIES ///////////////
		CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0));

		if(args->multiPhase){
#if !CUDA
			peridicBoundaries(n, m, r_f, b_f);
#else
			gpuBcPeriodic2D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, r_f_d, b_f_d,bcCount, cg_dir_d);
#endif

		} else{
			gpuBcInlet2D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					u0_d, v0_d, bcCount);
			gpuBcWall2D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					fColl_d, qCollapsed_d, bcCount);
			gpuBcOutlet2D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					u0_d, v0_d, bcCount);
		}

		//				if(args->multiPhase){
		//					CHECK(cudaMemcpy(r_f, r_f_d, SIZEFLT(m*n *9), cudaMemcpyDeviceToHost));
		//					CHECK(cudaMemcpy(b_f, b_f_d, SIZEFLT(m*n *9), cudaMemcpyDeviceToHost));
		//				}

		CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_BNDC] += cudatime;

		// UPDATE VELOCITY AND DENSITY
		CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0));
		if(args->multiPhase){
#if !CUDA
			updateMacroMP(n,m,u,v,r_rho, b_rho, r_f, b_f, rho, args->control_param,args->r_alpha, args->b_alpha,
					args->bubble_radius,st_error, iter,st_predicted);
#else
			gpuUpdateMacro2DCG<<<bpg1, tpb>>>(fluid_d, rho_d, u_d, v_d, r_f_d, b_f_d, r_rho_d, b_rho_d, p_in_d, p_out_d, num_in_d, num_out_d);

			//			CHECK(cudaMemcpy(r_rho, r_rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
			//			CHECK(cudaMemcpy(b_rho, b_rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
			//			updateSurfaceTension(r_rho,b_rho,args->control_param, st_predicted, st_error, iter,args->r_alpha, args->b_alpha, args->bubble_radius, n ,m);
			//gpu reduction is faster than serial surface tension
			p_in_mean = gpu_sum_h(p_in_d, p_in_d, ms) / gpu_sum_int_h(num_in_d, num_in_d, ms);
			p_out_mean = gpu_sum_h(p_out_d, p_out_d, ms) / gpu_sum_int_h(num_out_d, num_out_d, ms);
			st_error[iter] = calculateSurfaceTension(p_in_mean, p_out_mean,args->r_alpha, args->b_alpha, args->bubble_radius, st_predicted);
#endif
		}
		else gpuUpdateMacro2D<<<bpg1, tpb>>>(fluid_d, rho_d, u_d, v_d, bcMask_d,
				drag_d, lift_d, coordX_d, coordY_d, f_d);

		tInstant2 = clock();
		CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_MACR] += cudatime;

		// COMPUTE RESIDUALS
		CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0));
		FLOAT_TYPE r;
		if(args->multiPhase){
#if !CUDA
			CHECK(cudaMemcpy(r_f_d, r_f, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
			CHECK(cudaMemcpy(b_f_d, b_f, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
			CHECK(cudaMemcpy(r_fColl_d, r_fColl, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
			CHECK(cudaMemcpy(b_fColl_d, b_fColl, SIZEFLT(m*n *9), cudaMemcpyHostToDevice));
#endif
			r = computeResidual2D(r_f_d, r_fColl_d, temp9a_d, temp9b_d, m,n);

		}else
			r = computeResidual2D(f_d, fColl_d, temp9a_d, temp9b_d, m,n);
		if (r != r) {
			fprintf(stderr, "\nDIVERGENCE!\n");

			writeResiduals(residualsFilename, norm, dragSum, liftSum, m * n,
					iter + 1);
			cudaEventDestroy(start);
			cudaEventDestroy(stop);

			freeAllHost(hostArrays, sizeof(hostArrays) / sizeof(hostArrays[0]));
			freeAllGpu(gpuArrays, sizeof(gpuArrays) / sizeof(gpuArrays[0]));

			return 1; // ERROR!
		}
		norm[iter] = r;
		if (args->boundaryId > 0) {
			dragSum[iter] = computeDragLift2D(bcMask_d, drag_d, tempA_d,
					tempB_d, m, n, args->boundaryId);
			liftSum[iter] = computeDragLift2D(bcMask_d, lift_d, tempA_d,
					tempB_d, m, n, args->boundaryId);
		}

		CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_RESI] += cudatime;
		printf("Iterating... %d/%d (%3.1f %%)\r", iter + 1, args->iterations,
				(FLOAT_TYPE) (iter + 1) * 100
				/ (FLOAT_TYPE) (args->iterations));

		iter++; // update loop variable
		////////////// Autosave ///////////////

		if (iter == (args->autosaveEvery * autosaveIt)) {
			autosaveIt++;
			if (iter > args->autosaveAfter) {
				printf("autosave\n\n");
				//////////// COPY VARIABLES TO HOST ////////////////
				CHECK(cudaMemcpy(u, u_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
				CHECK(cudaMemcpy(v, v_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
				CHECK(cudaMemcpy(rho, rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
				if(args->multiPhase){
					CHECK(cudaMemcpy(r_rho, r_rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
					CHECK(cudaMemcpy(b_rho, b_rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
				}

				switch (args->outputFormat) {
				case CSV:
					sprintf(finalFilename, "%sFinalData.csv", inFn->result);
					break;
				case TECPLOT:
					sprintf(finalFilename, "%sFinalData.dat", inFn->result);
					break;
				case PARAVIEW:
					sprintf(finalFilename, "%sFinalData.vti", inFn->result);
					break;
				}

				tInstant1 = clock(); // Start measuring time
				WriteResults3D(finalFilename, nodeType,nodeX, nodeY, nodeZ, u, v, w, rho,
						nodeType, n, m, 1, args->outputFormat);
				tInstant2 = clock();
				taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1)
																																																																																																																										/ CLOCKS_PER_SEC;
			}
		}
	}     ////////////// END OF MAIN WHILE CYCLE! ///////////////
	tIterEnd = clock(); // End measuring time of main loop
	taskTime[T_ITER] = (FLOAT_TYPE) (tIterEnd - tIterStart) / CLOCKS_PER_SEC;

	clock_t tEnd = clock();
	taskTime[T_OALL] = (FLOAT_TYPE) (tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time
	taskTime[T_COLL] /= 1000;
	taskTime[T_STRM] /= 1000;
	taskTime[T_BNDC] /= 1000;
	taskTime[T_MACR] /= 1000;
	taskTime[T_RESI] /= 1000;
	fclose(logFile);
	writeEndLog(logFilename, taskTime);
	writeTimerLog(timeFilename, taskTime);
	writeResiduals(residualsFilename, norm, dragSum, liftSum, m * n,
			args->iterations);

#if CUDA
	//WRITE VARIABLES TO HOST
	CHECK(cudaMemcpy(u, u_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(v, v_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(rho, rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
	if(args->multiPhase){
		CHECK(cudaMemcpy(r_rho, r_rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(b_rho, b_rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
	}
#endif
	switch (args->outputFormat) {
	case CSV:
		sprintf(finalFilename, "%sFinalData.csv", inFn->result);
		break;
	case TECPLOT:
		sprintf(finalFilename, "%sFinalData.dat", inFn->result);
		break;
	case PARAVIEW:
		sprintf(finalFilename, "%sFinalData.vti", inFn->result);
		break;
	}
	if(args->multiPhase){
		WriteResultsMultiPhase(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,r_rho,b_rho, nodeType,
				n, m, 1, args->outputFormat);
		printf("Suface tension error: "FLOAT_FORMAT"\n", st_error[iter-1]);
		WriteArray("surface tension",st_error, args->iterations,1);
	}
	else
		WriteResults3D(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho, nodeType,
				n, m, 1, args->outputFormat);

	// Write information for user
	printf("\n\nLog was written to %s\n", logFilename);
	printf("Last autosave result can be found at %s\n", autosaveFilename);
	printf("residuals were written to %s\n", residualsFilename);
	printf("Profiling results were written to %s\n", timeFilename);
	printf("Final results were written to %s\n", finalFilename);

	//	compareTestFiles("./TestValues/CUDA/rpert.txt", "./TestValues/CUDA/rpert_gpu.txt");
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	freeAllHost(hostArrays, sizeof(hostArrays) / sizeof(hostArrays[0]));
	freeAllGpu(gpuArrays, sizeof(gpuArrays) / sizeof(gpuArrays[0]));
	return 0;
}
