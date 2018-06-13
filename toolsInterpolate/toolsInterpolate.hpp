//
//  toolsInterpolate.hpp
//  
//
//  Created by Darius Abramavicius on 10/5/12.
//  Copyright (c) 2012 VU. All rights reserved.
//
#pragma once

#ifndef _toolsInterpolate_hpp
#define _toolsInterpolate_hpp

class toolsInterpolate
{
    
public:
    
    toolsInterpolate();
    ~toolsInterpolate();

    // return value is the error code
    // 0 - everything is correct
    // 1 - some error core has been generated
    
    // simple linear interpolation
    int linint(           
           double *inX, 
           double *inY, 
           int n, 
           double x,  
           double& ret, 
           double& err);

    
    // simple polynomial interpolation
    int polint(
           double *inX, 
           double *inY, 
           int n, 
           double x,  
           double& ret, 
           double& err);


    // simple rational function interpolation
    int ratint(
           double* xa, 
           double* ya, 
           int n, 
           double x, 
           double& y, 
           double& dy);
    
    

    // Barycentric interpolation
    // this default setup uses d=n-1 interpolation
    int barint(
               double *inX, 
               double *inY, 
               int n, 
               double x,  
               double& ret, 
               double& err);
    // more advanced usage:
    // use initially barinit() , then barretval(), finally interclean()
    int barinit(double *inX, int n, int d);
    // this function calculates the weights
    // for the d order of Barycentric interpolation
    int barretval(double *inX, double *inY, double x, double& y);
    // this function returns interpolated value
    
    
private:
    double* tdata;  // some dataset depending on the algorithm
    int tint;   // some int for temporary storage
    
    void interclean();

};

#endif
