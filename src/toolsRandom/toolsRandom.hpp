#pragma once



#include"../complexv/complexv.hpp"
#include"../dvector3d/dvector3d.hpp"


class toolsRandom
{
    
public:
    
    // constructors
    toolsRandom(int); // with seed initialization
    toolsRandom();

        double RandLD(double k);//returns linear random from 0 to <k

        double RandED(double k);//returns exponential random

        double RandGD(double v, double s);//returns gaussian random

        int GaussianTrajectory(complexv* gtraj, int nump, double correlationL);

        // returns a randomly oriented unit length 3D vector
        // orientation is distributed uniformly on a sphere
    	dvector3d RandV3D()
    	{
    		// generating random vector:
    		double length;
    		double x,y,z;
    		do
    		{
    			x = -1.0 + RandLD(2.0);
    			y = -1.0 + RandLD(2.0);
    			z = -1.0 + RandLD(2.0);
    			length = (x*x + y*y + z*z);
    		}while ( length > 1.0 );
    		length = sqrt(length);
    		x = x / length;
    		y = y / length;
    		z = z / length;

    		return dvector3d(x,y,z);
    	}

};


