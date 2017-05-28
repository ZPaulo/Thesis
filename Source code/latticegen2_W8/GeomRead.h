#include <string>
#include <sstream>
#include <vector>
using namespace std;
#ifndef GEOMREAD_H
#define GEOMREAD_H


class GeomRead
{
    public:
        GeomRead();
        const std::vector<vector<double> > & GetNodes() const { return Nodes;}
        const std::vector<vector<double> > & GetConnectors() const { return Connectors;}
        int getNumNodes(){return NumNodes;};
        int getNumConnectors(){return NumConnectors;};
        static string filename;
    protected:
    private:
        vector<vector<double> > Nodes;
        vector<vector<double> > Connectors;
        int i;
        int j;
		int k;
        int NumNodes;
        int NumConnectors;
        int ID;
        double CoordX;
        double CoordY;
        double CoordZ;
        int Int1;
        int Point1, Point2, Point3, Point4;
        int Int2;
        int BCID;
        string BCtype;
};

#endif // GEOMREAD_H
