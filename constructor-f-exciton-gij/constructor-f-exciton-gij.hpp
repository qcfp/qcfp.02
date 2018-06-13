#pragma once



#include "../storage/storage.hpp"
#include "../numericalSD/numericalSD.hpp"


// class to process and manipulate Frenkel exciton Hamiltonian
// creates amplitudes to the correlation function in eigenstate basis
class constructor_f_exciton_gij
{

    
public:

	// constructors
	constructor_f_exciton_gij(storage<double>* ievecs);
        constructor_f_exciton_gij(storage<double>* ievecs,storage<double>* ievec2);
	constructor_f_exciton_gij(storage<double>& ievecs);
        constructor_f_exciton_gij(storage<double>& ievecs,storage<double>& ievec2);
	~constructor_f_exciton_gij();				// Destructor

public:

    
    // for uncorrelated fluctuations
    int GetGijAmplitudes11(storage<double>*, double&);
    int GetGijAmplitudes11(storage<double>*, double*);

    int GetRedAmplitudes11(storage<double>*, double&);
    int GetRedAmplitudes11(storage<double>*, double*);
    

    
    
    int GetGijAmplitudes21(storage<double>*, double&);
    int GetGijAmplitudes22(storage<double>*, double&);
    int GetGijAmplitudes21(storage<double>*, double*);
    int GetGijAmplitudes22(storage<double>*, double*);

    int GetRedAmplitudes22(storage<double>*, double&);
    int GetRedAmplitudes22(storage<double>*, double*);
    
    
    // for correlated fluctuations
    int GetGijAmplitudesC11(storage<double>& ret, double** cst);
    int GetRedAmplitudesC11(storage<double>& ret, double** cst);    
    int GetGijAmplitudesC21(storage<double>& ret, double** cst);
    int GetGijAmplitudesC22(storage<double>& ret, double** cst);
    int GetRedAmplitudesC22(storage<double>& ret, double** cst);
    
    // next return full tetradic matrices
    int GetFullRedAmplitudes11(storage<double>&, storage<double>&);
    int GetFullRedAmplitudes21(storage<double>&, storage<double>&);
    int GetFullRedAmplitudes22(storage<double>&, storage<double>&);

    // for uncorrelated fluctuations
    int GetFullRedAmplitudes11(double**** amps, double* coupls);
    int GetFullRedAmplitudes21(double**** amps, double* coupls);
    int GetFullRedAmplitudes22(double**** amps, double* coupls);

    // for correlated fluctuations
    int GetFullRedAmplitudes11(double**** amps, double** coupls);
    int GetFullRedAmplitudes21(double**** amps, double** coupls);
    int GetFullRedAmplitudes22(double**** amps, double** coupls);


    int GetFullAmplitudesMultimode(storage<double>& results, storage<double>& source);

private:

    storage<double>* evecs;
    storage<double>* evec2;
        

};


