#include <string>
#include <sstream>
#include <vector>
using namespace std;
#include "Vector3D.h"

Vector3D::Vector3D()
{

}


void Vector3D::arr(int z,int y, int x){
		array3D.resize(z);
		for (int i = 0; i < z; ++i)
        {
            array3D[i].resize(y);
			for (int j=0;j<y;j++){
				array3D[i][j].resize(x);
			}
        }

        for(int i=0;i<z;i++)
			{
			for(int j=0;j<y;j++)
				for (int k=0;k<x;k++)
				{
					array3D[i][j][k]=0;
				}	
			}
    }
