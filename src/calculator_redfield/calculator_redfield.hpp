#pragma once
#include"../dvector3d/dvector3d.hpp"
#include"../complexv/complexv.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"


// This class calculates Redfield relaxation matrices:
// Markovian secular and non-secular (full)
// Non-markovian full

// propagation is not done here
// use propagators instead
class calculator_redfield
{
public:
    
    // general flag
    int ready;

    // system parameters
    storage<double> evals;
    int numG;
    int numE;
    int numF;

    
    // bath:
    // temperature
    double tempr;

    // bath oscillators
    // storage< asymptoticLF<double> > spdens; // for a set of spectral densities 

    storage<asymptoticLF_complexv> mfun; // Mfunctions
    storage<asymptoticLF_complexv> gfun; // Gfunctions
    storage<asymptoticLF_complexv> gdun; // G 1-st derivative functions
    storage<asymptoticLF_complexv> cfun; // Correlation functions
    storage<double> kij; // amplitudes for all pairs of interstate couplings
    storage<double> gij; // amplitudes for all pairs of energy fluctuations
    storage<double> Mijkl; // full matrix of correlation amplitudes
    
    // transport rates:
    storage<double> rates; // Markovian secular population rates
    storage<complexv> dephasings; // Markovian secular dephasing rates
    storage<complexv> supermatrix; // Markovian full nonsecular rates
    storage<complexv> kernel; // nonMarkovian time dependent relaxation kernel
    storage<double> reorganizations;
    
    void Setup(); // this sets up from scratch
    
    
public:

    calculator_redfield();
    calculator_redfield(storage<double>& ielevs);
    calculator_redfield(storage<double>& ielevs, int numg, int nume, int numf);  
    void AddFluctuations(storage<asymptoticLF_complexv>& imfun,
                         storage<double>& ikij,
                         storage<double>& igij);
    void AddFluctuations(
         storage<asymptoticLF_complexv>& imfun,
         storage<asymptoticLF_complexv>& icfun,
         storage<asymptoticLF_complexv>& igfun,
         storage<asymptoticLF_complexv>& igdun,
         storage<double>& ikij,
         storage<double>& igij,
         storage<double>& iMijkl );
    //void AddMijkl(storage<complexv>& iM);
    void AddMijkl(storage<double>& iM);
    

    storage<double> AddTransportRates();
    storage<complexv> AddLifetimeDephasings();
    storage<complexv> AddPureDephasings();
    storage<complexv> GetRelaxationSuperoperatorM(int block);
    storage<complexv> GetMemoryKernel(double time, double& deltat, int& numT, int block); // this creates the full Born memory kernel
    storage<complexv> GetMemoryKern(int block); // this returns only block-organized M-function
    storage<double> GetReorganizations();
    storage<complexv> GetRelaxationSuperoperatorLindblad1(int block);
    storage<double> GetTransportRatesModRed();

    int flagLindblad;
    //int nonsecular;
    
};
