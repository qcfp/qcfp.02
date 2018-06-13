#include"toolsInterpolate2d.hpp"
#include"../constants/constants.hpp"


#include<iostream>
using std::cout;

toolsInterpolate2d::toolsInterpolate2d()
{
    // do nothing
}

toolsInterpolate2d::~toolsInterpolate2d()
{
    // do nothing
}


// this function creates the point value on a plain
// defined by the three points: v0 at (0,0),  v1 at (0,1),  v2 at (1,0)
//    *
//    **
//    ***
//    ****
double toolsInterpolate2d::GetTriangle(double& xarg, double& yarg,
        double& v0, double& v1, double& v2)
{
    return v0+(v2-v0)*xarg+(v1-v0)*yarg;
}



// this function creates the point value on a square
//    ****
//    ****
//    ****
//    ****
// for this purpose it divides the squate into four triangles
// and uses the triangles
double toolsInterpolate2d::GetSquare0011(double coox,double cooy,double* val)
{
//        this program interpolates a 2D point \"strength\" value in 3D:\n";
//        val[0] \"strength\" value at (0,0) coordinate\n";
//        val[1] \"strength\" value at (0,1) coordinate\n";
//        val[2] \"strength\" value at (1,1) coordinate\n";
//        val[3] \"strength\" value at (1,0) coordinate\n";
    
    // OK. So now assume that I have four points
    
    // now central value
    double valc=0.25*(val[0]+val[1]+val[2]+val[3]);
    
    // return value
    double retval = 0.0;
    
    // now I have four triangles
    
    // special coefficient cos(45)
    constants cst;
    const double cs45 = cst.cos45;
    
    
    //I check if I am inside the square?
    if(coox<0 || coox > 1 || cooy<0 || cooy>1){
        cout<<"Error: toolsInterpolate2d: the point is out of range\n";
        return 0.0;
    }
    
    // I need to determine in which triangle this point is located?
    // for that purpose I have two lines:
    // y=x
    // y=1-x
    if(cooy > 1-coox && cooy > coox)
    {
        // I am in upper triangle
        //cout<<"Upper\n";
     
        // new coordinates in rotated and expanded frame
        double nx =  (coox-0.5)+(cooy-0.5) ;
        double ny = -(coox-0.5)+(cooy-0.5) ;
        
        retval = GetTriangle(nx, ny, valc, val[1], val[2]);
    }
    else if(cooy > 1-coox && cooy < coox)
    {
        // I am in right triangle
        //cout<<"Right\n";

        // new coordinates in rotated and expanded frame
        double nx =  (coox-0.5)-(cooy-0.5) ;
        double ny =  (coox-0.5)+(cooy-0.5) ;

        retval = GetTriangle(nx, ny, valc, val[2], val[3]);
    }
    else if(cooy <= 1-coox && cooy < coox)
    {
        // I am in bottom triangle
        //cout<<"Bottom\n";

        // new coordinates in rotated and expanded frame
        double nx = -(coox-0.5)-(cooy-0.5) ;
        double ny =  (coox-0.5)-(cooy-0.5) ;

        retval = GetTriangle(nx, ny, valc, val[3], val[0]);
    }
    else 
    {
        // I am in left triangle
        //cout<<"Left\n";

        // new coordinates in rotated and expanded frame
        double nx = -(coox-0.5)+(cooy-0.5) ;
        double ny = -(coox-0.5)-(cooy-0.5) ;

        retval = GetTriangle(nx, ny, valc, val[0], val[1]);
    }

    return retval;
}

