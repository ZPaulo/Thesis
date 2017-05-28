// Mesh generation

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
using namespace std;
#include "GeomRead.h"
#include "Vector2D.h"
#include "Vector3D.h"
#include "MeshGen.h"
#include "math.h"
int MeshGen::n = 0;
int MeshGen::m = 0;
int MeshGen::h = 0;

MeshGen::MeshGen() {
	////////////////////// NUMBER OF NODES IN X-DIRECTION /////////////////////
	//n=512;//880;

	GeomNodes = Geo.GetNodes();

	Col1.resize(Geo.getNumNodes());
	Col2.resize(Geo.getNumNodes());
	Col3.resize(Geo.getNumNodes());

	for (int v = 0; v < Geo.getNumNodes(); v++) {
		Col1[v] = GeomNodes[v][0];
		Col2[v] = GeomNodes[v][1];
		Col3[v] = GeomNodes[v][2];

	}

	// find maximum and minimum x and y in order to mesh the entire geometry
	XMax = GetMax(Col1);
	YMax = GetMax(Col2);
	ZMax = GetMax(Col3);
	XMin = GetMin(Col1);
	YMin = GetMin(Col2);
	ZMin = GetMin(Col3);

	Xdomain = XMax - XMin;
	Ydomain = YMax - YMin;
	Zdomain = ZMax - ZMin;

	// find delta x and delta y
	Delta = Xdomain / (n - 1);
	m = ceil((Ydomain / Delta)) + 1;
	h = ceil((Zdomain / Delta)) + 1;

	// store all the nodes in a list "MeshNodes"
	Vector2D Mnod;
	Mnod.arr((n + 2) * (m + 2) * (h + 2), 5);
	MeshNodes = Mnod.getarray2D();
	IDnode = 0;

	// store all nodes ID in a matrix "IDmatrix"
	Vector3D IDm;
	IDm.arr((h + 2), (m + 2), (n + 2));

	IDmatrix = IDm.getarray3D();

	for (int k = 0; k < h + 2; k++) {
		for (int j = 0; j < m + 2; j++) {
			for (int i = 0; i < n + 2; i++) {

				MeshNodes[i + j * (n + 2) + k * (m + 2) * (n + 2)][0] = IDnode;
				MeshNodes[i + j * (n + 2) + k * (m + 2) * (n + 2)][1] = XMin
						+ Delta * (i - 1); // outer nodes are at least 0.5*delta from the fluid nodes
				MeshNodes[i + j * (n + 2) + k * (m + 2) * (n + 2)][2] = YMin
						+ Delta * (j - 1);
				MeshNodes[i + j * (n + 2) + k * (m + 2) * (n + 2)][3] = ZMin
						+ Delta * (k - 1);
				IDmatrix[k][j][i] = IDnode;
				IDnode++;
			}
		}
	}

}

double MeshGen::GetMax(vector<double> V) {

	double tempMax = V[0];
	int D = V.size();
	for (int i = 0; i < D; i++) {

		if (V[i] > tempMax)
			tempMax = V[i];
	}
	return tempMax;
}

double MeshGen::GetMin(vector<double> W) {

	double tempMin = W[0];
	int D = W.size();
	for (int i = 0; i < D; i++) {
		if (W[i] < tempMin)
			tempMin = W[i];
	}
	return tempMin;
}
