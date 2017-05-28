#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>       /* fabs */
#include "Vector2D.h"
#include "GeomRead.h"
#include "MeshGen.h"
using namespace std;
#include "Intersect.h"

Intersect::Intersect() {
	conn_maxX = 0;
	conn_maxY = 0;
	conn_maxZ = 0;
	conn_minX = 0;
	conn_minY = 0;
	conn_minZ = 0;
	ID = -1;
	maxX = 0;
	maxY = 0;
	maxZ = 0;
	minX = 0;
	minY = 0;
	minZ = 0;
	SlopeXY = 0;
	SlopeXZ = 0;
	SlopeYZ = 0;
	Stop = -1;
	Xinter = 0;
	Yinter = 0;
	Zinter = 0;
	IDmatrix = Mesh.GetIDmatrix();
	MeshNodes = Mesh.GetMeshNodes();
	Nodes = Geom.GetNodes();
	Connectors = Geom.GetConnectors();
	BC = -1;

}

bool Intersect::Inter(int IDNode1, int IDNode2) {

	Node1[0] = MeshNodes[IDNode1][1];
	Node2[0] = MeshNodes[IDNode2][1];
	Node1[1] = MeshNodes[IDNode1][2];
	Node2[1] = MeshNodes[IDNode2][2];
	Node1[2] = MeshNodes[IDNode1][3];
	Node2[2] = MeshNodes[IDNode2][3];

	maxX = max(Node1[0], Node2[0]);
	minX = min(Node1[0], Node2[0]);
	maxY = max(Node1[1], Node2[1]);
	minY = min(Node1[1], Node2[1]);
	maxZ = max(Node1[2], Node2[2]);
	minZ = min(Node1[2], Node2[2]);

	SlopeXY = (Node2[1] - Node1[1]) / (Node2[0] - Node1[0]);

	Stop = false;

	for (int i = 0; i < Geom.getNumConnectors(); i++) {
		//ONGRID HANDLING:

		//the same Y and X coords (dir 5, 6)
		if (Node1[1] == Node2[1] && Node1[0] == Node2[0]) {

			if (Connectors[i][6] != 3 && Connectors[i][6] != 2) { // select non YZ and XZ faces

				Intersection[0] = Node1[0];
				Intersection[1] = Node1[1];

				if (Connectors[i][6] == 1) { //select XY faces

					Intersection[2] = Nodes[Connectors[i][2] - 1][2];

					if ((Intersection[2] <= maxZ) && (Intersection[2] >= minZ)
							&& (Intersection[0] <= Connectors[i][10])
							&& (Intersection[0] >= Connectors[i][11])
							&& (Intersection[1] <= Connectors[i][12])
							&& (Intersection[1] >= Connectors[i][13])) {
						if ((Intersection[0] != Node2[0] || Node1[0] == Node2[0])
								&& (Intersection[1] != Node2[1]
										|| Node1[1] == Node2[1])
								&& (Intersection[2] != Node2[2]
										|| Node1[2] == Node2[2])) {
							Stop = true;
							ID = Connectors[i][0];
							BC = Connectors[i][5];
							Xinter = Intersection[0];
							Yinter = Intersection[1];
							Zinter = Intersection[2];
							break;
						}

					}
				} else {
					/* T0DO SLOPING FACES
					 //select non vertical and non horizontal connectors
					 Intersection[0]=Nodes[Connectors[i][2]-1][0]+(Node1[1]-Nodes[Connectors[i][2]-1][1])/Connectors[i][5];
					 if ((Intersection[0]<=maxX)&&(Intersection[0]>=minX)&&(Intersection[1]<=Connectors[i][8])&&(Intersection[1]>=Connectors[i][9])){
					 Stop=true;
					 ID=Connectors[i][0];
					 BC=Connectors[i][3];
					 Xinter=Intersection[0];
					 Yinter=Intersection[1];
					 }
					 */
				}
			}

		}

		//the same Y and Z coords (dir 1, 2)
		if (Node1[1] == Node2[1] && Node1[2] == Node2[2]) {

			if (Connectors[i][6] != 1 && Connectors[i][6] != 2) { // select non XZ and XY faces

				Intersection[2] = Node1[2];
				Intersection[1] = Node1[1];

				if (Connectors[i][6] == 3) { //select XY faces

					Intersection[0] = Nodes[Connectors[i][2] - 1][0];

					if ((Intersection[0] <= maxX) && (Intersection[0] >= minX)
							&& (Intersection[1] <= Connectors[i][12])
							&& (Intersection[1] >= Connectors[i][13])
							&& (Intersection[2] <= Connectors[i][14])
							&& (Intersection[2] >= Connectors[i][15])) {
						if ((Intersection[0] != Node2[0] || Node1[0] == Node2[0])
								&& (Intersection[1] != Node2[1]
										|| Node1[1] == Node2[1])
								&& (Intersection[2] != Node2[2]
										|| Node1[2] == Node2[2])) {
							Stop = true;
							ID = Connectors[i][0];
							BC = Connectors[i][5];
							Xinter = Intersection[0];
							Yinter = Intersection[1];
							Zinter = Intersection[2];
							break;
						}
					}
				} else {
					/*
					 T0DO SLOPING FACES
					 //select non vertical and non horizontal connectors
					 Intersection[0]=Nodes[Connectors[i][2]-1][0]+(Node1[1]-Nodes[Connectors[i][2]-1][1])/Connectors[i][5];
					 if ((Intersection[0]<=maxX)&&(Intersection[0]>=minX)&&(Intersection[1]<=Connectors[i][8])&&(Intersection[1]>=Connectors[i][9])){
					 Stop=true;
					 ID=Connectors[i][0];
					 BC=Connectors[i][3];
					 Xinter=Intersection[0];
					 Yinter=Intersection[1];
					 }
					 */
				}
			}

		}

		//the same X and Z coords (dir 3, 4)
		if (Node1[0] == Node2[0] && Node1[2] == Node2[2]) {

			if (Connectors[i][6] != 1 && Connectors[i][6] != 3) { // select non XZ and XY faces

				Intersection[2] = Node1[2];
				Intersection[0] = Node1[0];

				if (Connectors[i][6] == 2) { //select XY faces

					Intersection[1] = Nodes[Connectors[i][2] - 1][1];

					if ((Intersection[1] <= maxY) && (Intersection[1] >= minY)
							&& (Intersection[0] <= Connectors[i][10])
							&& (Intersection[0] >= Connectors[i][11])
							&& (Intersection[2] <= Connectors[i][14])
							&& (Intersection[2] >= Connectors[i][15])) {
						if ((Intersection[0] != Node2[0] || Node1[0] == Node2[0])
								&& (Intersection[1] != Node2[1]
										|| Node1[1] == Node2[1])
								&& (Intersection[2] != Node2[2]
										|| Node1[2] == Node2[2])) {
							Stop = true;
							ID = Connectors[i][0];
							BC = Connectors[i][5];
							Xinter = Intersection[0];
							Yinter = Intersection[1];
							Zinter = Intersection[2];
							break;
						}
					}
				} else {
					/* T0DO SLOPING FACES
					 //select non vertical and non horizontal connectors
					 Intersection[0]=Nodes[Connectors[i][2]-1][0]+(Node1[1]-Nodes[Connectors[i][2]-1][1])/Connectors[i][5];
					 if ((Intersection[0]<=maxX)&&(Intersection[0]>=minX)&&(Intersection[1]<=Connectors[i][8])&&(Intersection[1]>=Connectors[i][9])){
					 Stop=true;
					 ID=Connectors[i][0];
					 BC=Connectors[i][3];
					 Xinter=Intersection[0];
					 Yinter=Intersection[1];
					 }
					 */
				}
			}

		}

		//if sloping (on XY plane) (dir 7,8,9,10)
		if ((Node1[0] != Node2[0]) && (Node1[1] != Node2[1])
				&& (Node1[2] == Node2[2])) {

			if (Connectors[i][6] != 1) { //select non XY Faces

				if (Connectors[i][6] == 2) { //select XZ faces
					//printf("here!\n");
					Intersection[2] = Node1[2];
					if (Node2[0] >= Node1[0] && Node2[0] < Connectors[i][10]) {
						Intersection[0] = Node1[0]
								+ fabs(
										Nodes[Connectors[i][2] - 1][1]
												- Node1[1]);

					} else {
						if (Node2[0] < Node1[0]
								&& Node2[0] >= Connectors[i][11]) {
							Intersection[0] = Node1[0]
									- fabs(
											Nodes[Connectors[i][2] - 1][1]
													- Node1[1]);
						}

					}
					if (Node1[1] <= Connectors[i][12]
							&& Node1[0] >= Connectors[i][11]
							&& Node1[1] == Connectors[i][12]) {
						Intersection[0] = Node1[0];
					}

					Intersection[1] = Nodes[Connectors[i][2] - 1][1];

					if ((Intersection[0] <= maxX) && (Intersection[0] >= minX)
							&& (Intersection[1] <= maxY)
							&& (Intersection[1] >= minY)
							&& (Intersection[2] <= Connectors[i][14])
							&& (Intersection[2] >= Connectors[i][15])) {
						if ((Intersection[0] != Node2[0] || Node1[0] == Node2[0])
								&& (Intersection[1] != Node2[1]
										|| Node1[1] == Node2[1])
								&& (Intersection[2] != Node2[2]
										|| Node1[2] == Node2[2])) {
							Stop = true;
							ID = Connectors[i][0];
							BC = Connectors[i][5];
							Xinter = Intersection[0];
							Yinter = Intersection[1];
							Zinter = Intersection[2];
							break;
						}
					}
				} else {
					if (Connectors[i][6] == 3) { //select YZ faces

						Intersection[2] = Node1[2];

						if (Node2[1] > Node1[1]
								&& Node2[1] <= Connectors[i][12]) {
							Intersection[1] = Node1[1]
									+ fabs(
											Nodes[Connectors[i][2] - 1][0]
													- Node1[0]);
						} else {
							if (Node2[1] < Node1[1]
									&& Node2[1] >= Connectors[i][13]) {
								Intersection[1] = Node1[1]
										- fabs(
												Nodes[Connectors[i][2] - 1][0]
														- Node1[0]);
							}
						}
						if (Node1[1] <= Connectors[i][12]
								&& Node1[1] >= Connectors[i][13]
								&& Node1[0] == Connectors[i][10]) {
							Intersection[1] = Node1[1];
						}
						Intersection[0] = Nodes[Connectors[i][2] - 1][0];

						if ((Intersection[0] <= maxX)
								&& (Intersection[0] >= minX)
								&& (Intersection[1] <= maxY)
								&& (Intersection[1] >= minY)
								&& (Intersection[2] <= Connectors[i][14])
								&& (Intersection[2] >= Connectors[i][15])) {
							if ((Intersection[0] != Node2[0]
									|| Node1[0] == Node2[0])
									&& (Intersection[1] != Node2[1]
											|| Node1[1] == Node2[1])
									&& (Intersection[2] != Node2[2]
											|| Node1[2] == Node2[2])) {
								Stop = true;
								ID = Connectors[i][0];
								BC = Connectors[i][5];
								Xinter = Intersection[0];
								Yinter = Intersection[1];
								Zinter = Intersection[2];
								break;
							}

						}
					}
				}
			}
		}

		//if sloping (on XZ plane) (dir 11,12,13,14)
		if ((Node1[0] != Node2[0]) && (Node1[2] != Node2[2])
				&& (Node1[1] == Node2[1])) {

			if (Connectors[i][6] != 2) { //select non XZ Faces

				if (Connectors[i][6] == 1) { //select XY faces

					Intersection[1] = Node1[1];
					if (Node2[0] > Node1[0] && Node2[0] <= Connectors[i][10]) {
						Intersection[0] = Node1[0]
								+ fabs(
										Nodes[Connectors[i][2] - 1][2]
												- Node1[2]);
					} else {
						if (Node2[0] < Node1[0]
								&& Node2[0] >= Connectors[i][11])
							Intersection[0] = Node1[0]
									- fabs(
											Nodes[Connectors[i][2] - 1][2]
													- Node1[2]);
					}
					if (Node1[0] <= Connectors[i][10]
							&& Node1[0] >= Connectors[i][11]
							&& Node1[2] == Connectors[i][14]) {
						Intersection[0] = Node1[0];
					}
					Intersection[2] = Nodes[Connectors[i][2] - 1][2];

					if ((Intersection[0] <= maxX) && (Intersection[0] >= minX)
							&& (Intersection[2] <= maxZ)
							&& (Intersection[2] >= minZ)
							&& (Intersection[1] <= Connectors[i][12])
							&& (Intersection[1] >= Connectors[i][13])) {
						if ((Intersection[0] != Node2[0] || Node1[0] == Node2[0])
								&& (Intersection[1] != Node2[1]
										|| Node1[1] == Node2[1])
								&& (Intersection[2] != Node2[2]
										|| Node1[2] == Node2[2])) {

							Stop = true;
							ID = Connectors[i][0];
							BC = Connectors[i][5];
							Xinter = Intersection[0];
							Yinter = Intersection[1];
							Zinter = Intersection[2];
							break;
						}

					}
				} else {
					if (Connectors[i][6] == 3) { //select YZ faces

						Intersection[1] = Node1[1];

						if (Node2[2] > Node1[2]
								&& Node2[2] <= Connectors[i][14]) {
							Intersection[2] = Node1[2]
									+ fabs(
											Nodes[Connectors[i][2] - 1][0]
													- Node1[0]);
						} else {
							if (Node2[2] < Node1[2]
									&& Node2[2] >= Connectors[i][15]) {
								Intersection[2] = Node1[2]
										- fabs(
												Nodes[Connectors[i][2] - 1][0]
														- Node1[0]);
							}
						}
						if (Node1[2] <= Connectors[i][14]
								&& Node1[2] >= Connectors[i][15]
								&& Node1[0] == Connectors[i][10]) {
							Intersection[2] = Node1[2];
						}
						Intersection[0] = Nodes[Connectors[i][2] - 1][0];

						if ((Intersection[0] <= maxX)
								&& (Intersection[0] >= minX)
								&& (Intersection[2] <= maxZ)
								&& (Intersection[2] >= minZ)
								&& (Intersection[1] <= Connectors[i][12])
								&& (Intersection[1] >= Connectors[i][13])) {
							if ((Intersection[0] != Node2[0]
									|| Node1[0] == Node2[0])
									&& (Intersection[1] != Node2[1]
											|| Node1[1] == Node2[1])
									&& (Intersection[2] != Node2[2]
											|| Node1[2] == Node2[2])) {
								Stop = true;
								ID = Connectors[i][0];
								BC = Connectors[i][5];
								Xinter = Intersection[0];
								Yinter = Intersection[1];
								Zinter = Intersection[2];
								break;
							}

						}
					}
				}
			}

		}
		//if sloping (on YZ plane) (dir 15,16,17,18)
		if ((Node1[2] != Node2[2]) && (Node1[1] != Node2[1])
				&& (Node1[0] == Node2[0])) {

			if (Connectors[i][6] != 3) { //select non YZ Faces

				if (Connectors[i][6] == 1) { //select XY faces

					Intersection[0] = Node1[0];
					if (Node2[1] > Node1[1] && Node2[1] <= Connectors[i][12]) {
						Intersection[1] = Node1[1]
								+ fabs(
										Nodes[Connectors[i][2] - 1][2]
												- Node1[2]);
					} else {
						if (Node2[1] < Node1[1]
								&& Node2[1] >= Connectors[i][13]) {
							Intersection[1] = Node1[1]
									- fabs(
											Nodes[Connectors[i][2] - 1][2]
													- Node1[2]);
						}
					}
					if (Node1[1] <= Connectors[i][12]
							&& Node1[1] >= Connectors[i][13]
							&& Node1[2] == Connectors[i][14]) {
						Intersection[1] = Node1[1];
					}
					Intersection[2] = Nodes[Connectors[i][2] - 1][2];

					if ((Intersection[2] <= maxZ) && (Intersection[2] >= minZ)
							&& (Intersection[1] <= maxY)
							&& (Intersection[1] >= minY)
							&& (Intersection[0] <= Connectors[i][10])
							&& (Intersection[0] >= Connectors[i][11])) {
						if ((Intersection[0] != Node2[0] || Node1[0] == Node2[0])
								&& (Intersection[1] != Node2[1]
										|| Node1[1] == Node2[1])
								&& (Intersection[2] != Node2[2]
										|| Node1[2] == Node2[2])) {
							Stop = true;
							ID = Connectors[i][0];
							BC = Connectors[i][5];
							Xinter = Intersection[0];
							Yinter = Intersection[1];
							Zinter = Intersection[2];
							break;
						}

					}
				} else {
					if (Connectors[i][6] == 2) { //select XZ faces

						Intersection[0] = Node1[0];

						if (Node2[2] > Node1[2]
								&& Node2[2] <= Connectors[i][14]) {
							Intersection[2] = Node1[2]
									+ fabs(
											Nodes[Connectors[i][2] - 1][1]
													- Node1[1]);
						} else {
							if (Node2[2] < Node1[2]
									&& Node2[2] >= Connectors[i][15]) {
								Intersection[2] = Node1[2]
										- fabs(
												Nodes[Connectors[i][2] - 1][1]
														- Node1[1]);
							}
						}
						if (Node1[2] <= Connectors[i][14]
								&& Node1[2] >= Connectors[i][15]
								&& Node1[1] == Connectors[i][12]) {
							Intersection[2] = Node1[2];
						}
						Intersection[1] = Nodes[Connectors[i][2] - 1][1];

						if ((Intersection[2] <= maxZ)
								&& (Intersection[2] >= minZ)
								&& (Intersection[1] <= maxY)
								&& (Intersection[1] >= minY)
								&& (Intersection[0] <= Connectors[i][10])
								&& (Intersection[0] >= Connectors[i][11])) {
							if ((Intersection[0] != Node2[0]
									|| Node1[0] == Node2[0])
									&& (Intersection[1] != Node2[1]
											|| Node1[1] == Node2[1])
									&& (Intersection[2] != Node2[2]
											|| Node1[2] == Node2[2])) {
								Stop = true;
								ID = Connectors[i][0];
								BC = Connectors[i][5];
								Xinter = Intersection[0];
								Yinter = Intersection[1];
								Zinter = Intersection[2];
								break;
							}

						}
					}
				}
			}

		}
	}

	return Stop;
}

