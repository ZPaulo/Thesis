#ifndef INTERSECT_H
#define INTERSECT_H


class Intersect
{
    public:

        Intersect();
        bool Inter(int IDNode1,int IDNode2);
        double getXinter(){return Xinter;};
        double getYinter(){return Yinter;};
		double getZinter(){return Zinter;};
        int getBC(){return BC;};
        int getID(){return ID;};
        //bool getn(){return Stop;};
    protected:
    private:

        double Node1[3];
        double Node2[3];
        double SlopeXY;
		double SlopeXZ;
		double SlopeYZ;
        MeshGen Mesh;
        GeomRead Geom;
        vector<vector<vector<double> > > IDmatrix;
        vector<vector<double> > Nodes;
        vector<vector<double> > Connectors;
        vector<vector<double> > MeshNodes;
        double Intersection[3];
        double maxY;
        double minY;
		double minX;
        double maxX;
        double minZ;
		double maxZ;
        double Yinter;
        double Xinter;
		double Zinter;
        double conn_maxY;
        double conn_minY;
        double conn_maxX;
        double conn_minX;
		double conn_maxZ;
        double conn_minZ;
        int BC;
        int ID;
        bool Stop;
};

#endif // INTERSECT_H
