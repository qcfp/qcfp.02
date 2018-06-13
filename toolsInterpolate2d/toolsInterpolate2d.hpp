//
//  toolsInterpolate2d.hpp
//  
//
#pragma once

class toolsInterpolate2d
{
    
public:
    
    toolsInterpolate2d();
    ~toolsInterpolate2d();


    double GetSquare0011(double coox,double cooy,double* val);
    
    double GetTriangle(double& xarg, double& yarg,
        double& v0, double& v1, double& v2);
    

};

