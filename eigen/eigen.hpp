#pragma once



#include "../storage/storage.hpp"
#include "../complexv/complexv.hpp"
//#include"../asymptoticLF/asymptoticLF.hpp"


// class to process and manipulate eigenvalue problem
class eigen
{

    
public:

	// constructors
	eigen(int);  // used for making eigenvectors in rows (int !=0 )
	eigen();
	~eigen();				// Destructor

public:

    int GetEigenSolution(storage<double>& imatrix, storage<double>& evals, storage<double>& evecs);
    int GetEigenSolution(storage<double>* imatrix, storage<double>* evals, storage<double>* evecs);
    int GetEigenSolution(storage<double>& imatrix, storage<complexv>& evals, storage<complexv>& evecsL,  storage<complexv>& evecsR);
    int GetEigenSolution(storage<complexv>& imatrix, storage<complexv>& evals, storage<complexv>& evecsL,  storage<complexv>& evecsR);
private:
    
    int ltype;
    //double smallepsilon;
    
    
    /////////////////////////////
    // this part uses the standard LAPACK routines
    

};


