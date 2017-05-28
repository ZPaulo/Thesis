// Starting from the seed point, the neighbors in the eight directions are checked using the method "Inter" in the class "INTERSECT"

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "Vector2D.h"
#include "Vector3D.h"
#include "GeomRead.h"
#include "MeshGen.h"
#include "Intersect.h"
#include "FluidTracking.h"
#include <iomanip>

std::string FluidTracking::filename = "";
double FluidTracking::XSeed = 0;
double FluidTracking::YSeed = 0;
double FluidTracking::ZSeed = 0;
int FluidTracking::GeometryComplexity = 1;

FluidTracking::FluidTracking() {
	ofstream BCconn(filename.c_str());

	//     D3Q19 LATTICE CONFIGURATION

	// 			  Z= 1 LAYER
	//                15       
	//                |     
	//                |   
	//                | 
	//        12 - - -5 - - - 11
	//                | 
	//                |  
	//                |     
	//                16   

	// 			  Z=0 LAYER
	//        8       3       7
	//          \     |     /
	//            \   |   /
	//              \ | /
	//        2 - - - 0 - - - 1
	//              / | \
	//            /   |   \
	//          /     |     \
	//        10       4      9

	// 			  Z= -1 LAYER
	//                17       
	//                |     
	//                |   
	//                | 
	//        14 - - -6 - - - 13
	//                | 
	//                |  
	//                |     
	//                18      

	// BCconnectors.dat
	// |    Xindex    |   Yindex	| Zindex	|	ID lattice   |   BC type |   X intersection   |   Y intersection   |	Z intersection   | BCID

	GeomRead GEOM;
	MeshGen MESH;
	Intersect INTER;
	int cx3D[19] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0,
			0 };
	int cy3D[19] = { 0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1,
			-1 };
	int cz3D[19] = { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1,
			-1 };
	IDmatrix = MESH.GetIDmatrix();
	MeshNodes = MESH.GetMeshNodes();
	Connectors = GEOM.GetConnectors();

	m = MESH.m;
	n = MESH.n;
	h = MESH.h;
	int X = -1, Y = -1, Z = -1;
	Vector3D a;
	a.arr((h + 2), (m + 2), (n + 2));
	Fluidmatrix = a.getarray3D();

	Percentage = m * n * h / 100;
	perc = 0;

	// SEED POINT/////
	// XSeed=0.5;//70;
	// YSeed=0.5;//27;
	////////////////////
	// identify the closest node
	Distance = sqrt(pow(XSeed, 2) + pow(YSeed, 2) + pow(ZSeed, 2));
	for (int i = 0; i < (m + 2) * (n + 2) * (h + 2); i++) {
		//cout<<i<<endl;
		if (sqrt(
				pow(MeshNodes[i][1] - XSeed, 2)
						+ pow(MeshNodes[i][2] - YSeed, 2)
						+ pow(MeshNodes[i][3] - ZSeed, 2)) < Distance) {
			Distance = sqrt(
					pow(MeshNodes[i][1] - XSeed, 2)
							+ pow(MeshNodes[i][2] - YSeed, 2)
							+ pow(MeshNodes[i][3] - ZSeed, 2));
			IDseed = MeshNodes[i][0];
		}
	}
	// append the seed node in the list of all the fluid nodes
	Fluid.push_back(IDseed);
	MeshNodes[IDseed][4] = 1;

	k = 0;

	// for each node listed in "Fluid" check the 18 directions

	while (k < Fluid.size()) {

		Zindex = (int) (Fluid[k] / ((n + 2) * (m + 2)));
		Yindex = (int) ((Fluid[k] - Zindex * (n + 2) * (m + 2)) / (n + 2));
		Xindex = (int) ((Fluid[k] - Zindex * (n + 2) * (m + 2)
				- Yindex * (n + 2)));

		// check if the segments linking the node with the surrounding nodes intersect a boundary.
		// if it does it the surrounding node is written in the connectors file, otherwise the surrounding node is appended in the "Fluid" list
		// the intersection is checked with the class "Intersect"
		if (Xindex > 0 && Xindex < n + 1 && Yindex > 0 && Yindex < m + 1
				&& Zindex > 0 && Zindex < h + 1) {

			for (int dir = 1; dir < 19; dir++) {
				if (INTER.Inter(IDmatrix[Zindex][Yindex][Xindex],
						IDmatrix[Zindex + cz3D[dir]][Yindex + cy3D[dir]][Xindex
								+ cx3D[dir]]) == false) {

					if (Fluidmatrix[Zindex + cz3D[dir]][Yindex + cy3D[dir]][Xindex
							+ cx3D[dir]] == 0) {
						Fluid.push_back(
								IDmatrix[Zindex + cz3D[dir]][Yindex + cy3D[dir]][Xindex
										+ cx3D[dir]]);
						MeshNodes[IDmatrix[Zindex + cz3D[dir]][Yindex
								+ cy3D[dir]][Xindex + cx3D[dir]]][4] = 1;
						Fluidmatrix[Zindex + cz3D[dir]][Yindex + cy3D[dir]][Xindex
								+ cx3D[dir]] = 1;

					}
				} else {
					if (GeometryComplexity == 1) {

						if (Fluidmatrix[Zindex + cz3D[dir]][Yindex + cy3D[dir]][Xindex
								+ cx3D[dir]] == 0) {
							WriteConn(BCconn, Xindex, Yindex, Zindex, dir,
									INTER.getBC(), INTER.getXinter(),
									INTER.getYinter(), INTER.getZinter(),
									INTER.getID());
							BCFluid.push_back(IDmatrix[Zindex][Yindex][Xindex]);
						}
					}

				}

			}

		}
		k++;
		if (k > perc * Percentage) {
			if (GeometryComplexity == 2) {
				cout << "\rComplete:  " << k * 100 / n / m / h / 2 << "%"
						<< flush;
			} else {
				cout << "\rComplete:  " << k * 100 / n / m / h << "%" << flush;
			}

			perc++;
		}

	}
	if (GeometryComplexity == 2) {
		k = 0;
		while (k < Fluid.size()) {

			Zindex = (int) (Fluid[k] / ((n + 2) * (m + 2)));
			Yindex = (int) ((Fluid[k] - Zindex * (n + 2) * (m + 2)) / (n + 2));
			Xindex = (int) ((Fluid[k] - Zindex * (n + 2) * (m + 2)
					- Yindex * (n + 2)));

			for (int dir = 1; dir < 19; dir++) {
				if (Xindex > 0 && Xindex < n + 1 && Yindex > 0 && Yindex < m + 1
						&& Zindex > 0 && Zindex < h + 1) {

					if (INTER.Inter(IDmatrix[Zindex][Yindex][Xindex],
							IDmatrix[Zindex + cz3D[dir]][Yindex + cy3D[dir]][Xindex
									+ cx3D[dir]]) != false) {

						if (Fluidmatrix[Zindex + cz3D[dir]][Yindex + cy3D[dir]][Xindex
								+ cx3D[dir]] == 0) {
							WriteConn(BCconn, Xindex, Yindex, Zindex, dir,
									INTER.getBC(), INTER.getXinter(),
									INTER.getYinter(), INTER.getZinter(),
									INTER.getID());
							BCFluid.push_back(IDmatrix[Zindex][Yindex][Xindex]);
						}
					}
				}
			}

			k++;
			if (k > (perc - 50) * Percentage) {
				cout << "\rComplete:  " << k * 100 / n / m / h / 2 + 50 << "%"
						<< flush;
				perc++;
			}
		}
	}

	cout << "\rComplete:  " << 100 << "%" << flush;
	cout << endl;

	ofstream Outfil("xtecplot.dat");

	Outfil << "title = \"sample mesh\"" << endl;
	Outfil << "variables = \"x\", \"y\", \"z\" " << endl;
	Outfil << "zone i=" << n + 2 << ", j=" << m + 2 << ", k=" << h + 2 << endl;
	for (int k = 0; k < (h + 2) * (m + 2) * (n + 2); k++) {
		Outfil << MeshNodes[k][1] << "    " << MeshNodes[k][2] << "    "
				<< MeshNodes[k][3] << endl;
	}

}

void FluidTracking::WriteConn(ofstream& fstream, int Xindex, int Yindex,
		int Zindex, int dir, int BC, float Xinter, float Yinter, float Zinter,
		int BCID) {

	fstream << fixed << setprecision(6) << Xindex - 1 << "\t" << Yindex - 1
			<< "\t" << Zindex - 1 << "\t" << dir << "\t" << BC << "\t" << Xinter
			<< "\t" << Yinter << "\t" << Zinter << "\t" << BCID << endl;
}

