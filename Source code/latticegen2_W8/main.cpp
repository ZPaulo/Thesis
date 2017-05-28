#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include "Vector2D.h"
#include "GeomRead.h"
#include "MeshGen.h"
#include "Intersect.h"
#include "FluidTracking.h"
#include <ctime>
using namespace std;



int main(int argc, char** argv)
{
    int m,n,h;
    clock_t begin=clock();

    if (argc != 8) {
           cout << "Usage: meshgen <input_filename> <NoOfNodes in X dir>  <xseed> <yseed> <zseed> <output filename> <1- simple, 2 - complex geometry>\n";
   		cout<<argc<<endl;
           return 1;
       }

    GeomRead::filename = argv[1];
    MeshGen::n = atoi(argv[2]);
    FluidTracking::XSeed = atof(argv[3]);
    FluidTracking::YSeed = atof(argv[4]);
	FluidTracking::ZSeed = atof(argv[5]);
	string OutputFilename = argv[6];
    FluidTracking::filename = OutputFilename+"_BC.dat";
    FluidTracking::GeometryComplexity=atoi(argv[7]);

//%%%%%%%%%%%%%% CALL GEOMETRY READING CONTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GeomRead GEOM;
	

//%%%%%%%%%%%%%% CALL MESH GENERATION CONTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MeshGen MESH;

//%%%%%%%%%%%%%% CALL FLUIDTRACKING CONTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FluidTracking Fluid;

	n=MESH.n;
    m=MESH.m;
	h=MESH.h;
  
//%%%%%%%%%%%%%% WRITE NODE FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//  i   |   j   |   X  |    Y   |   Fluid/Solid |
	string temp=OutputFilename+"_Node.dat";
    ofstream D2Nod(temp.c_str());

	for(int k=1;k<h+1;k++){
		for(int j=1;j<m+1;j++){
			for(int i=1;i<n+1;i++){

				D2Nod<<i-1<<"\t"<<j-1<<"\t"<<k-1<<fixed<<setprecision(6)<<"\t"<<Fluid.GetMeshNodes()[i+j*(n+2)+k*(m+2)*(n+2)][1]<<"\t"<<Fluid.GetMeshNodes()[i+j*(n+2)+k*(m+2)*(n+2)][2]<<"\t"<<Fluid.GetMeshNodes()[i+j*(n+2)+k*(m+2)*(n+2)][3]<<"\t"<<(int)Fluid.GetMeshNodes()[i+j*(n+2)+k*(m+2)*(n+2)][4]<<endl;
			}
		}
	}
	 clock_t end = clock();
	  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	  cout<<"work took: "<<elapsed_secs<<"s"<<endl;
    return 0;
}
