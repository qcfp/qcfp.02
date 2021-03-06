/*
 * File:   propagatorMemory.hpp
 * Author: dariusa
 *
 * Created on October 29, 2015, 10:09 AM
 */

#ifndef PROPAGATORMEMORY_HPP
#define	PROPAGATORMEMORY_HPP


// ***********************
// master equation with memory
//
// the translationary invariant time-depenant kernel is prepared by the user
// the class makes
// propagation
// usage :
// Init, Convolute, Update,  and Exit



#pragma once


#include"../storage/storage.hpp"
#include"../complexv/complexv.hpp"
#include"../propagatorM/propagatorM.hpp"
#include"../interpolationF/interpolationF.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"

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
    void Initialize(storage<complexv>& ikernel,
                    storage<double>& iomegas_memory,
										storage<double>& iomegas_reorganizations_memory,
                    storage<complexv>& iDMM,
                    storage<double>& itimesinternal);
    void Initialize(propagatorMemory& iprop);
    void InitializeHistory(propagatorMemory& iprop);

    // this function calculates current variables when derivs is given
    void Update( storage<complexv>& der, double timestep);

    // this function calculates derivatives when kernel is given
    void Convolute(storage<complexv>& der);
		void ConvoluteTC2(storage<complexv>& der);
		void ConvoluteTCL2(storage<complexv>& der);

    // this function makes the next step when derivatives are given
    storage<complexv> Propagate(storage<double>& times);


    // accessory functions to convert to/from Liouville space
    storage<complexv> Convert3DTo5D(storage<complexv> input,int num4,int num3,int num2,int num1);
    storage<complexv> Convert5DTo3D(storage<complexv> input);
    storage<complexv> Convert2DTo3D(storage<complexv> input,int num2,int num1);
    storage<complexv> Convert3DTo2D(storage<complexv> input);
    storage<complexv> Convert2DTo4D(storage<complexv> input,int num4,int num3,int num2,int num1);
    storage<complexv> Convert4DTo2D(storage<complexv> input);
    storage<complexv> Convert1DTo2D(storage<complexv> input,int num2,int num1);
    storage<complexv> Convert2DTo1D(storage<complexv> input);

    // absolute history times
    // notice: time goes in reverse direction
    storage<double> timesinternal;
    double internaltimeS;
    int internaltimeN;

		int approachTC;
		// this is the method o. propagation either
		// time convolution 1
		// or local time convolutionless 2



    // block 0:
    // kernel
    storage<complexv>* kernel;

    // additional history tuner when not time translationary
    storage<double> omegas_memory;
		storage<double> omegas_reorganizations_memory;

    //history matrix 1
    storage<complexv>* DMM;


    // notice that blocks are not used with cfun.... not clear if they make any use

    // block 1:
    storage<complexv>* kernel1;
    storage<double> omegas_memory1;
    storage<complexv>* DMM1;

    // block 2:
    storage<complexv>* kernel2;
    storage<double> omegas_memory2;
    storage<complexv>* DMM2;

    // entry points of different blocks
    double manifold0End;
    double manifold1End;
    double manifold2End;
};


#endif	/* PROPAGATORMEMORY_HPP */
