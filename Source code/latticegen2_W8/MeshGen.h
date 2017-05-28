#include <string>
#include <sstream>
#include <vector>
#include <fstream>
using namespace std;
#include "GeomRead.h"
#include "Vector2D.h"
#ifndef MESHGEN_H
#define MESHGEN_H


class MeshGen
{
    public:
        MeshGen();
        const std::vector<vector<double> > & GetMeshNodes() const { return MeshNodes;}
        const std::vector<vector<vector<double> > > & GetIDmatrix() const { return IDmatrix;}
        double getDelta(){return Delta;};
        double getXmax(){return XMax;};
        double getYmax(){return YMax;};
        double getZmax(){return ZMax;};
        double getXmin(){return XMin;};
        double getYmin(){return YMin;};
        double getZmin(){return ZMin;};

        static int n; //number of nodes in x direction
		static int m; //number of nodes in y direction
		static int h; //number of nodes in z direction
    protected:
    private:
        GeomRead Geo;
        vector<vector<double> > MeshNodes;
        vector<vector<double> > GeomNodes;
        vector<vector<vector<double> > > IDmatrix;
        vector<double> Col1;
        vector<double> Col2;
		 vector<double> Col3;
        double XMax;
        double YMax;
		double ZMax;
        double XMin;
        double YMin;
		double ZMin;
        double Xdomain;
        double Ydomain;
		double Zdomain;
        double Delta;
        int IDnode;
        double GetMax(vector<double> V);
        double GetMin(vector<double> W);
};

#endif // MESHGEN_H
