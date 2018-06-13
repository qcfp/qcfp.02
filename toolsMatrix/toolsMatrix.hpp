#pragma once
// forward declarations to avoid confusing file inclusion problems
template <typename T> class storage;
#include"../complexv/complexv.hpp"
#include"../toolsRandom/toolsRandom.hpp"


class toolsMatrix
{
public:

	//default constructor;
	toolsMatrix();
	toolsMatrix(int);
	toolsMatrix(toolsRandom* itrd);
    ~toolsMatrix();
    

    
    
    // setting random vectors
    void VecRandLD(storage<double>&, double k);//sets elements linear random from 0 to <k
    void VecRandED(storage<double>&, double k);//sets elements exponential random
    void VecRandGD(storage<double>&, double v, double s);//sets elements gaussian random

    
    
    
    
    
    
	void Set0(storage<double>&);
	void Set1(storage<double>&);
    void AddVecToDiag(storage<double>& matr,storage<double>& vec);


	
    
    
    
    
    
    
    
    
    //operators
	// assignment
	void Prod (storage<double>&,storage<double>&,storage<double>&);
	void Add (storage<double>&,storage<double>&,storage<double>&);
	void InverseC (storage<complexv>&,storage<complexv>&);
        void InverseC (complexv* ret, complexv* mat, int n);
    
    
    
    
    // members
    toolsRandom* trd;
    int trdallocated;
    

};

