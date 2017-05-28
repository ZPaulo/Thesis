#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#ifndef VECTOR3D_H
#define VECTOR3D_H


class Vector3D
{
    public:
        Vector3D();
        void arr(int X,int Y, int Z);
        const std::vector<vector<vector<double> > > & getarray3D() const { return array3D;}
    protected:
    private:
        vector<vector<vector<double> > > array3D;
};

#endif // 3DVECTOR_H
