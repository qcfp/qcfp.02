#pragma once



#include "../storage/storage.hpp"
#include "../dvector3d/dvector3d.hpp"
#include "../toolsRandom/toolsRandom.hpp"
//#include"../asymptoticLF/asymptoticLF.hpp"


// class to disorder properties
class constructor_disorder
{

    
public:

	// constructors
	constructor_disorder( int iset);
	constructor_disorder();
   
	~constructor_disorder();// Destructor

public:

    // add Gaussian intra-diagonal disorder
    // onto the matrix
    int DiagIntra(storage<double>& matr, double sigma); 
    
    // add Gaussian inter-diagonal disorder
    // onto the matrix
    int DiagInter(storage<double>& matr, double sigma); 
    
private:
    
    toolsRandom* processor;
    
    int sval;
};


