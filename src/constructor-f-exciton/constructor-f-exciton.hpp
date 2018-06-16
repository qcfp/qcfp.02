#pragma once



#include "../storage/storage.hpp"
#include "../dvector3d/dvector3d.hpp"
#include "../dtensor3x3/dtensor3x3.hpp"
//#include"../asymptoticLF/asymptoticLF.hpp"

// when doing eigenvectors
// we use 2D  arrays 
// evecs[a][b]
// here a is the site number
// b is the eigenstate number


// class to process and manipulate Frenkel exciton Hamiltonian
class constructor_f_exciton
{

    
public:

	// constructors
	constructor_f_exciton(
                          storage<double>& ham,
                          storage<dvector3d>& dips
                          );
    constructor_f_exciton(
                          int lev,
                          storage<double>& ham,
                          storage<dvector3d>& dips
                          );
    constructor_f_exciton(storage<double>& ham);
    constructor_f_exciton(int lev, storage<double>& ham);
	constructor_f_exciton();
	~constructor_f_exciton();				// Destructor

public:

    int GetEvals(storage<double>&);
    
    // GetEvecs returns [site][estate]
    int GetEvecs(storage<double>&);

    // GetEvecsES returns [estate][site]
    int GetEvecsES(storage<double>&);
    int GetEdips(storage<dvector3d>&);
    int GetEvals2(storage<double>&);
    int GetEvecs2(storage<double>&);
    int GetStateEvecs2();
    int GetEdips2(storage<dvector3d>&);
    
    int GetHam2(storage<double>&);
	int GetHamiltonian2LocalSorters(storage<int>& sl,storage<int>& sr);
    int GetDips(storage<dvector3d>&);
    int GetDips2(storage<dvector3d>&);
    
    int AddKCouplings(storage<double>& iK);
    int AddADipoles(storage<dvector3d>& id);

	int AddDipolePositions(storage<dvector3d>& ip);
	int AddMagneticDipoles(storage<dvector3d>& im);
	int GetETens(storage<dtensor3x3>& ten);
	int GetEMags(storage<dvector3d>& mag);
	int GetETens2(storage<dtensor3x3>& ten);
	int GetEMags2(storage<dvector3d>& mag);

	int GetTens(storage<dtensor3x3>& ten);
	int GetMags(storage<dvector3d>& mag);
	int GetTens2(storage<dtensor3x3>& ten);
	int GetMags2(storage<dvector3d>& mag);

    
    void reformatF12();
    void reformatF21();
    
private:
    
    storage<double> evals;
    storage<double> evecs;    
    storage<dvector3d> edips;
    
    storage<double> evals2;
    storage<double> evecs2;
    storage<dvector3d> edips2;
    
    
    storage<double> hamloc;
    storage<dvector3d> dipsloc;

    storage<double> hamloc2;
    storage<dvector3d> dipsloc2;

	// representer of the hamloc2 indexer into the pairs of sites
    	storage<int> manifold2sorter_u;
	storage<int> manifold2sorter_l;

    storage<double> hamK;
    storage<dvector3d> dips2; // anharmonic dipoles 

    storage<dvector3d>  lpos; // local coordinates of transitions
    storage<dvector3d>  lmag; // local magnetic dipoles
    storage<dtensor3x3>  lten; // local quadrupole tensors of transitions
    storage<dtensor3x3>  lten2; // eigentensors of quadrupole transitions
    storage<dvector3d>  lmag2; // eigente magnetic transition dipoles

    storage<dtensor3x3>  eten; // eigentensors of quadrupole transitions
    storage<dvector3d>  emag; // eigente magnetic transition dipoles
    storage<dtensor3x3>  eten2; // eigentensors of quadrupole transitions
    storage<dvector3d>  emag2; // eigente magnetic transition dipoles

    
    int levels; // how many levels in a site (2 or 3)

    
private:
    int Init();
    int MakeEigens1();
    int MakeEigens2();
    int MakeHamLoc2();
    int MakeEdips1();
    int MakeEdips2();
    int Makedipsloc2();
    int Makedipsloc();
};


