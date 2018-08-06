#pragma once
#include"../dvector3d/dvector3d.hpp"
#include"../dtensor3x3/dtensor3x3.hpp"
#include"../complexv/complexv.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
#include"../interpolationF2d/interpolationF2d.hpp"
#include"../toolsIO/toolsIO.hpp"

#include<string>
#include<fstream>

using std::string;
using std::ifstream;
using std::ofstream;

// no allocation / deallocation occurs here
// all pointers are passed from outside
class communicator_3rd
{
public:

	// all information - text input
	stringstream inputfile;
    
    // general flag
    int ready;
    
    string complexity; // for the system
    string bath_complexity;
    string experiment_complexity;

    int band; 
    // 1 - one exciton properties
    // 2 - 1 and 2 exciton properties 

    int numT; // total number of levels
    
    //////////
    // level system:
    // manifolds
    bool eigensys;
    int numG;
    int numE;
    int numF;
    int nstates;
    // system parameters
    storage<double> evals;
    storage<dvector3d> edips;
    storage<complexv> dephasings;
    storage<double> grPops; // ground manifold populations
    storage<double> reorganizations;

	// more advanced eigensystem parameters:
	storage<dtensor3x3> eten; // quadrupole transition tensors (Angstroms * dipoles)
	storage<dvector3d> ed_m; // eigenstate magnetic moments of the transitions (??? units will be cleared in the future)
    
    
    /////////
    // excitonic system:
    bool bosonic;
    storage<double> ham;
    storage<double> anh;
    storage<dvector3d> dip;
    storage<dvector3d> dip2Corrections;
    storage<double> coupl; // amplitudes for all couplings [osc][chromophore]
    int excitonictransformed;
    storage<double> evec0; // for consistency
    storage<double> evec1;
    storage<double> evec2;

	// more advanced excitonic parameters:
	storage<dvector3d> pos; // positions of excitonic transitions (Angstroms)
	storage<dvector3d> d_m; // magnetic moments of the transitions (??? units will be cleared in the future)
    	storage<int> ham2sorter_l;
	storage<int> ham2sorter_r;

    
    // bath:
    // bath oscillators
    double tempr; // temperature
    storage<asymptoticLF<double> > spdf; // spectral density
    storage<asymptoticLF_complexv> gfun; // gfunctions
    storage<asymptoticLF_complexv> mfun; // Mfunctions
    storage<asymptoticLF_complexv> cfun; // Cfunctions
    storage<double> gij; // amplitudes for all pairs of energy levels
    storage<double> kij; // amplitudes for all pairs of interstate couplings
    storage<double> zijkl; // all fluctuating correlation coefficients
    // transport rates:
    storage<double> ratesg;
    storage<double> ratese;
    double transportstep;
    int transport;
    int nonsecular;
    
    double propagationparameterStep;
    double propagationparameterNump;

    // Experiments
    // interaction configurations:
    int coherentGSBK1;
    int transportGSBK1;
    int coherentESEK1;
    int transportESEK1;
    int coherentESAK1;
    int transportESAK1;
    int coherentGSBK2;
    int transportGSBK2;
    int coherentESEK2;
    int transportESEK2;
    int coherentESAK2;
    int transportESAK2;
    int coherentES1K3;
    int coherentES2K3;
    
    
        // experiment parameters
    int nump;
    double ifre1;
    double ffre1;
    double ifre2;
    double ffre2;
    double ifre3;
    double ffre3;
    double tf1;
    double tf2;
    double tf3;
    double ti1;
    double ti2;
    double ti3;
    
    dvector3d vecE1;
    dvector3d vecE2;
    dvector3d vecE3;
    dvector3d vecE4;

	// more advanced parameters of the excitations (wavevectors - only directions)
	dvector3d veck1;
	dvector3d veck2;
	dvector3d veck3;
	dvector3d veck4;

	// more advanced parameters of the excitations (excitation central frequencies - in wavenumbers cm-1)
	double xomega1;
	double xomega2;
	double xomega3;
	double xomega4;

	// more advanced parameters of the excitations (excitation bandwidths - sigmas of Gaussians in wavenumbers cm-1)
	double xsigma1;
	double xsigma2;
	double xsigma3;
	double xsigma4;

    int averagingtype;
    
    double naturallinewidth;

    // signal data
    int simType;
    string outputstring;
    string interactionpattern;

    interpolationF2d<complexv> signal;
    interpolationF<complexv> signal1d;
    
    void Setup(); // this sets up from scratch
    void Delete(); // this sets up from scratch
    
    
    communicator_3rd* GetThis()
    {
        return this;
    }
    
    
public:
    
    communicator_3rd();
    
   void SetEFields(dvector3d& ivE1,
                   dvector3d& ivE2,
                   dvector3d& ivE3,
                   dvector3d& ivE4);
    void SetEFields(dvector3d& ivE1,
                   dvector3d& ivE2);
   
    void SetAveraging(int& avera);
    void SetTemperature(double&);
    void SetupManifolds(int, int);
    void SetupManifolds(int, int, int);
    void AddTransport();
	void AddTransport(storage<double>& iratesG,storage<double>& iratesE,double tstep);
	void AddDephasings();
	void AddDephasingsAll();

    void FinalizeESSystem(int nonmarkovian);
    void SwitchConfigurations(
                              int icoherentGSBK1,
                              int itransportGSBK1,
                              int icoherentESEK1,
                              int itransportESEK1,
                              int icoherentESAK1,
                              int itransportESAK1,
                              int icoherentGSBK2,
                              int itransportGSBK2,
                              int icoherentESEK2,
                              int itransportESEK2,
                              int icoherentESAK2,
                              int itransportESAK2,
                              int icoherentES1K3,
                              int icoherentES2K3  );
    
    void AddPattern(string istr);


    void SwitchKI();
    void SwitchKII();
    void SwitchKIII();
    void AddFluctuations(storage<asymptoticLF_complexv> mfun,storage<double> ikij,storage<double> igij);
    void AddFluctuations(storage<asymptoticLF_complexv> gfun,storage<asymptoticLF_complexv> mfun,storage<double> ikij,storage<double> igij);
    
	void MakeESSystem();
	void MakeESSystem(int flag);
    
    virtual void ReadSpecific(ifstream& ifs);
//    virtual void o_ReadSpecific(Open& input);


    void ReadSystem(ifstream& f);
    void ReadBath(ifstream& f);
    void ReadExperiment(ifstream& f);
    void ReadExperimentLin(ifstream& f);
    void ReadExperiment3rd(ifstream& f);
    virtual void ReadFile(ifstream& f);
    virtual void ReadFile(string& f);
    void Publish(ofstream& cco);
    void Publish(string& cco);
    void Publish2DRe(ofstream& ofs);
    void Publish2DIm(ofstream& ofs);

    void def_ReadFile(string& filename);
//    virtual void o_ReadFile(string& filename);
    
    virtual void PublishSpecial4(ofstream& ofs){};


	//int LookUpAndReadSquareMatrix(string keyword,string label, int numl, int numr, storage<double>& arr);
	//int LookUpAndReadTriangularMatrix(string keyword,string label, int numl, storage<double>& arr);

	int LookUpAndReadBathSpectraldensities(string keyword,string label, int numBosc);

};
