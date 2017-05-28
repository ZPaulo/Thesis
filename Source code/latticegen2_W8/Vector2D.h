#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#ifndef VECTOR2D_H
#define VECTOR2D_H


class Vector2D
{
    public:
        Vector2D();
        void arr(int HEIGHT,int WIDTH);
        //vector<vector<double> > array2D;
        const std::vector<vector<double> > & getarray2D() const { return array2D;}
        //const vector<double>& get_array2D() const {return array2D;}
    protected:
    private:
        vector<vector<double> > array2D;
};

#endif // 2DVECTOR_H
