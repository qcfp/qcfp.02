#pragma once
#include"../complexv/complexv.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
#include"../interaction/interaction.hpp"

class feinman2sideddiagram1
{
private:
    
    // no allocation is performed in this object: all parameters are copies
    // no deallocation is done
    
    // energy (imaginary part is for the dephasing)
    complexv freq;
    
    // amplitudes of transitions (dipoles)
    dvector3d *ds;
    
    // vectors of fields (dipoles)
    dvector3d *es;
    
    // amplitude of the diagram
    complexv damplitude;
    
    // bath:
    asymptoticLF_complexv* gij;
    double* gijAmplitude;
    double* gijAmplitude2;
    double* gijAmplitude3;
    double* gijAmplitude4;
    int numOsc;
    
public:
    feinman2sideddiagram1();
    feinman2sideddiagram1(complexv& ifreq, interaction& itr);
    feinman2sideddiagram1(complexv& ifreq, dvector3d* ids, dvector3d* ies,int& averaging);
    feinman2sideddiagram1(complexv& ifreq, double& amplitude);
    ~feinman2sideddiagram1();
    
    void assignLineshape(asymptoticLF_complexv* gij,double* gijAmplitude, int& num);
    void assignLineshape(asymptoticLF_complexv* gij,double* gijA1,double* gijA2,double* gijA3,double* gijA4, int& num);
    void propagate(interpolationF<complexv>&);
    void propagateC(interpolationF<complexv>&);
    complexv calculate(double&);
    //void printout();
    
};

