#include <string>
#include <sstream>
#include <vector>
using namespace std;
#include "Vector2D.h"

Vector2D::Vector2D()
{

}


void Vector2D::arr(int HEIGHT,int WIDTH){
       array2D.resize(HEIGHT);
        for (int i = 0; i < HEIGHT; ++i)
        {
            array2D[i].resize(WIDTH);
        }

        for(int k=0;k<HEIGHT;k++)
        {
        for(int j=0;j<WIDTH;j++)
            array2D[k][j]=0;
        }

    }
