#pragma once

// this project interfaces with FFTW3 standard GNU library
#include<fftw3.h>
#include"../interpolationF/interpolationF.hpp"
#include"../interpolationF2d/interpolationF2d.hpp"
#include"../complexv/complexv.hpp"


class toolsFFT
{

public: 
    toolsFFT();
	~toolsFFT();
    
    
    // one dimensional Fourier transforms
	void executeP(complexv* result,complexv* source,int);
	void executeN(complexv* result,complexv* source,int);
	interpolationF<complexv> executeP(interpolationF<complexv>& source);
	interpolationF<complexv> executeN(interpolationF<complexv>& source);
    
    
    // two dimensional Fourier transforms
	void executeP(complexv** result,complexv** source,int,int);
	void executeN(complexv** result,complexv** source,int,int);
	interpolationF2d<complexv> executeP(interpolationF2d<complexv>& source);
	interpolationF2d<complexv> executeN(interpolationF2d<complexv>& source);
        
        
        void SwapSides(interpolationF<complexv>& asymresFFTf);
        void SwapSides(interpolationF2d<complexv>& asymresFFTf);
    
    
    void cleanFFTSpace();


private:
    int nump;
    int nump2;
 int locked;
 int dimension;
 void executeGen(int& selection);
    
    
 fftw_complex* FFTdatI;
 fftw_complex* FFTdatO;
 fftw_plan fplanP;
 fftw_plan fplanN;

 void prepareFFTSpace(int np);
void prepareFFTSpace(int d2,int d1);


};
