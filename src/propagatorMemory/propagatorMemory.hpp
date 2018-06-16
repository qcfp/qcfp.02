/* 
 * File:   propagatorMemory.hpp
 * Author: dariusa
 *
 * Created on October 29, 2015, 10:09 AM
 */

#ifndef PROPAGATORMEMORY_HPP
#define	PROPAGATORMEMORY_HPP


#pragma once


#include"../storage/storage.hpp"
#include"../complexv/complexv.hpp"
#include"../propagatorM/propagatorM.hpp"
#include"../interpolationF/interpolationF.hpp"

class propagatorMemory:
public propagatorM
{
public:
    
    // simply designed to solve equation 
    // d p / dt = - int_0^{\infty} ds kernel(s) * p(t-s)
    
	propagatorMemory();
    
	~propagatorMemory();
    // usage : create the kernel (kernel) and the initial vector with history (DMM)
    // Propagate

    
    // this function initializes the algorithm
    void Initialize(storage<complexv>& kernel, 
                    storage<double>& omegas_memory, 
                    storage<complexv>& DMM, 
                    storage<double>& timesinternal);
    void Initialize(propagatorMemory& prop);
    void InitializeHistory(propagatorMemory& prop);
    
    // this function calculates updates current variables when derivs is given
    void Update( storage<complexv>& der,  // current derivative values
                double timestep);
    
    // this function calculates derivatives when kernel is given
    void Convolute(storage<complexv>& der
                    );

    storage<complexv> Propagate(storage<double>& times);
    

    // accessory functions
    storage<complexv> Convert3DTo5D(storage<complexv> input,int num4,int num3,int num2,int num1);
    storage<complexv> Convert5DTo3D(storage<complexv> input);
    storage<complexv> Convert2DTo3D(storage<complexv> input,int num2,int num1);
    storage<complexv> Convert3DTo2D(storage<complexv> input);
    storage<complexv> Convert2DTo4D(storage<complexv> input,int num4,int num3,int num2,int num1);
    storage<complexv> Convert4DTo2D(storage<complexv> input);
    storage<complexv> Convert1DTo2D(storage<complexv> input,int num2,int num1);
    storage<complexv> Convert2DTo1D(storage<complexv> input);


    storage<double> timesinternal;


    storage<complexv>* kernel;
    storage<double> omegas_memory;
    storage<complexv>* DMM;
    
    storage<complexv>* kernel1;
    storage<double> omegas_memory1;
    storage<complexv>* DMM1;

    storage<complexv>* kernel2;
    storage<double> omegas_memory2;
    storage<complexv>* DMM2;
    
    double manifold0End;
    double manifold1End;
    double manifold2End;
    
};


#endif	/* PROPAGATORMEMORY_HPP */

