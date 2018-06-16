#pragma once
#include"../dvector3d/dvector3d.hpp"
#include"../complexv/complexv.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
#include"../communicator_3rd/communicator_3rd.hpp"


// no allocation / deallocation occurs here
// all pointers are passed from outside
class calculator_abs_secular_cumulant :
public communicator_3rd
{
public:
    

//    void AddLineshapes(storage< asymptoticLF_complexv>& igfun,
//                       storage< asymptoticLF_complexv>& imfun,
//                       storage<double>& igij,storage<double>& ikij, double itempr);

    calculator_abs_secular_cumulant();
    calculator_abs_secular_cumulant(string& ifname);
//    calculator_abs_secular_cumulant(
//        storage<double>& ielevs,
//        storage<dvector3d>& iedips,
//        storage<complexv>& idephasings);
    ~calculator_abs_secular_cumulant();


//    // this function uses all input as pointers not allocating anything
//    interpolationF<complexv> Launch(double ifre, double ffre, int nump);
    
    // this function allocated everything
    virtual void Launch();

    
//    void ReadSpecific(ifstream& ifs);
    int liouville_pathway[3];
    int noComplexLifetimes;
	int markovian;

};

