#ifndef FLUIDTRACKING_H
#define FLUIDTRACKING_H


class FluidTracking
{
    public:
        FluidTracking();
        const std::vector<vector<double> > & GetMeshNodes() const { return MeshNodes;}
        static double XSeed;
        static double YSeed;
		static double ZSeed;
        static std::string filename;
        static int GeometryComplexity;
    protected:
    private:
        void WriteConn(ofstream& fstream, int Xindex, int Yindex, int Zindex, int dir, int BC, float Xinter, float Yinter, float Zinter, int BCID);
        vector<vector<vector<double> > > IDmatrix;
        vector<vector<vector<double> > > Fluidmatrix;
        vector<vector<double> > MeshNodes;
        vector<vector<double> > Connectors;
        vector<int> Fluid;
        vector<int> BCFluid;
        double Distance;
        double IDseed;
        int Xindex;
        int Yindex;
		int Zindex;
        int k;
        int m;
        int n;
		int h;
        int w;
        int perc;
        int Percentage;
};

#endif // FLUIDTRACKING_H
