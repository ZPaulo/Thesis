//
// STAR-CD files reading

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
using namespace std;
#include "GeomRead.h"
#include "Vector2D.h"

string GeomRead::filename = "";

GeomRead::GeomRead() {
	CoordX = -1;
	CoordY = -1;
	CoordZ = -1;
	ID = -1;
	Int1 = -1;
	Int2 = -1;
	BCID = -1;
	Point1 = -1;
	Point2 = -1;
	Point3 = -1;
	Point4 = -1;
	i = 0;
	k = 0;
	j = 0;
//%%%%%%%%%%%%%%%%%%% OPEN STAR CD FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	string nodeFile = "STARCDfiles/" + filename + ".vrt";
	string connFile = "STARCDfiles/" + filename + ".bnd";
	ifstream CountNodes(nodeFile.c_str());
	ifstream GeometryNodes(nodeFile.c_str());
	ifstream CountConnectors(connFile.c_str());
	ifstream GeometryConnectors(connFile.c_str());
	if (!CountNodes.is_open() || !CountConnectors.is_open()) {
		printf("Wrong input file name!\n");

	}
	while (CountNodes >> ID >> CoordX >> CoordY >> CoordZ) {
		i++;
	}

//store node info in a matrix "Nodes"
	NumNodes = i;
	Vector2D Nod;
	Nod.arr(i, 3);
	Nodes = Nod.getarray2D();
	i = 0;

	while (GeometryNodes >> ID >> CoordX >> CoordY >> CoordZ) {
		Nodes[i][0] = CoordX;
		Nodes[i][1] = CoordY;
		Nodes[i][2] = CoordZ;
		i++;

	}

	while (CountConnectors >> Int1 >> Point1 >> Point2 >> Point3 >> Point4
			>> BCID >> Int2 >> BCtype) {
		j++;
	}

	//store connectors info in a matrix "Connectors"
	NumConnectors = j;
	Vector2D Con;
	Con.arr(j, 16);
	Connectors = Con.getarray2D();

	j = 0;
	while (GeometryConnectors >> Int1 >> Point1 >> Point2 >> Point3 >> Point4
			>> BCID >> Int2 >> BCtype) {

		// set ID of the BC (given by pointwise)
		Connectors[j][0] = BCID;
		// set ID of the first point
		Connectors[j][1] = Point1;
		// set ID of the second point
		Connectors[j][2] = Point2;
		// set ID of the third point
		Connectors[j][3] = Point3;
		// set ID of the fourth point
		Connectors[j][4] = Point4;
		// set BC type
		if (BCtype == "WALL") {
			Connectors[j][5] = 1;
		}
		if (BCtype == "INLE") {
			Connectors[j][5] = 2;
		}
		if (BCtype == "OUTL") {
			Connectors[j][5] = 3;
		}
		if (BCtype == "CYCL") {
			Connectors[j][5] = 4;
		}
		if (BCtype == "SYMP") {
			Connectors[j][5] = 5;
		}

		/* 2d CASE
		 // set segment slope
		 if (Nodes[Point1-1][0]==Nodes[Point2-1][0]){
		 Connectors[j][6]=1;//vertical connector
		 }
		 else{

		 if (Nodes[Point1-1][1]==Nodes[Point2-1][1]){
		 Connectors[j][6]=-1;//horizontal connector
		 }
		 else{
		 Connectors[j][6]=0; // neither orizontal or vertical
		 Connectors[j][7]=(Nodes[Point2-1][1]-Nodes[Point1-1][1])/(Nodes[Point2-1][0]-Nodes[Point1-1][0]); // slope
		 }

		 }
		 */
		if (Nodes[Point1 - 1][0] == Nodes[Point2 - 1][0]
				&& Nodes[Point2 - 1][0] == Nodes[Point3 - 1][0]
				&& Nodes[Point3 - 1][0] == Nodes[Point4 - 1][0]) {
			Connectors[j][6] = 3; //YZ plane
		} else {

			if (Nodes[Point1 - 1][1] == Nodes[Point2 - 1][1]
					&& Nodes[Point2 - 1][1] == Nodes[Point3 - 1][1]
					&& Nodes[Point3 - 1][1] == Nodes[Point4 - 1][1]) {
				Connectors[j][6] = 2; //XZ plane
				//printf("%f\n",Connectors[j][1]);
			} else {
				if (Nodes[Point1 - 1][2] == Nodes[Point2 - 1][2]
						&& Nodes[Point2 - 1][2] == Nodes[Point3 - 1][2]
						&& Nodes[Point3 - 1][2] == Nodes[Point4 - 1][2]) {
					Connectors[j][6] = 1;                //XY plane
				} else {

					Connectors[j][6] = 0; // neither orizontal or vertical
					Connectors[j][7] = 0; // reserve T0DO
					Connectors[j][8] = 0; // reserve T0DO
					Connectors[j][9] = 0; // reserve T0DO
				}
			}
		}

		Connectors[j][10] = max(Nodes[Connectors[j][2] - 1][0],
				Nodes[Connectors[j][1] - 1][0]); //max x
		Connectors[j][10] = max(Connectors[j][10],
				Nodes[Connectors[j][3] - 1][0]); //max x
		Connectors[j][10] = max(Connectors[j][10],
				Nodes[Connectors[j][4] - 1][0]); //max x

		Connectors[j][11] = min(Nodes[Connectors[j][2] - 1][0],
				Nodes[Connectors[j][1] - 1][0]); //min x
		Connectors[j][11] = min(Connectors[j][11],
				Nodes[Connectors[j][3] - 1][0]); //min x
		Connectors[j][11] = min(Connectors[j][11],
				Nodes[Connectors[j][4] - 1][0]); //min x

		Connectors[j][12] = max(Nodes[Connectors[j][2] - 1][1],
				Nodes[Connectors[j][1] - 1][1]); //max y
		Connectors[j][12] = max(Connectors[j][12],
				Nodes[Connectors[j][3] - 1][1]); //max y
		Connectors[j][12] = max(Connectors[j][12],
				Nodes[Connectors[j][4] - 1][1]); //max y

		Connectors[j][13] = min(Nodes[Connectors[j][2] - 1][1],
				Nodes[Connectors[j][1] - 1][1]); //min y
		Connectors[j][13] = min(Connectors[j][13],
				Nodes[Connectors[j][3] - 1][1]); //min y
		Connectors[j][13] = min(Connectors[j][13],
				Nodes[Connectors[j][4] - 1][1]); //min y

		Connectors[j][14] = max(Nodes[Connectors[j][2] - 1][2],
				Nodes[Connectors[j][1] - 1][2]); //max z
		Connectors[j][14] = max(Connectors[j][14],
				Nodes[Connectors[j][3] - 1][2]); //max z
		Connectors[j][14] = max(Connectors[j][14],
				Nodes[Connectors[j][4] - 1][2]); //max z

		Connectors[j][15] = min(Nodes[Connectors[j][2] - 1][2],
				Nodes[Connectors[j][1] - 1][2]); //min z
		Connectors[j][15] = min(Connectors[j][15],
				Nodes[Connectors[j][3] - 1][2]); //min z
		Connectors[j][15] = min(Connectors[j][15],
				Nodes[Connectors[j][4] - 1][2]); //min z

		j++;
	}

}
