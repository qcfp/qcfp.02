#pragma once



#include"../complexv/complexv.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"


// class to process and manipulate spectral-density data
// at one fixed temperature
class numericalSD
{

	// the spectral densities are specified as 
	// as a functions of energy [in reciprocal centimeters]

public:

	double globalgamma;

	// constructors
    numericalSD();
	numericalSD(double dtemperature);
	numericalSD(double dtemperature,asymptoticLF<double>&);

	~numericalSD();				// Destructor

public:
    
    // constants that can be updated depending on the units
    // that are being used
    // consider exp(omega*time)
    // it is replaced by 
    // exp(omega*time/const5305)
    double const5305;
    // so all frequencies are used in "omega" units
    // and all times are used in "time" units
    // default is 1
    
    // the temperature wherever used is used in "omega" units.
    double const0_695; // this is Boltzmann constant
    // it can be redefined so that whatever units can be used
    // default is 1
    
    // the constant
    double const3000;
    // is tied with others.
    // it must be 2*pi*const5305
    
    void updateUnits(string istring);
    // updates unit system
    // "default"
    // "cm-fs-K"
    // are possible
    
	void ReadSD(string);
        void SetSD(double* idat,double dx,int nump);
        void SetSD(asymptoticLF<double>&);


	double GetSD(double);
        asymptoticLF<double> GetSD();// returns the whole function

	complexv GetGf(double);
	asymptoticLF_complexv GetGf();// returns the whole function
	complexv GetGfD(double, int); // up to second derivative
        complexv GetGfA(double);// gets g function from Acc functions

	// returns lineshape function first derivative at infinite time (when all correlations have decayed). 
	complexv GetGfD1Inf();

    complexv GetMf(double omega);
	complexv GetMfp(double omega);
	complexv GetMfn(double omega);
    asymptoticLF_complexv GetMfn();// returns the whole function
    asymptoticLF_complexv GetMf();// returns the whole function

	// accessory functions
	complexv GetA1(double);
	complexv GetA2(double);
	complexv GetA3(double);
        
	complexv GetCf(double);
	asymptoticLF_complexv GetCf();
    asymptoticLF_complexv GetGfD1f();




	void printMBAR(string fname);
	void printACCESSORY(string fname);
	void printGFUN(string fname);
	void printGDERIV1(string fname);
	void printGDERIV2(string fname);
	void printCORELF(string fname);

	void resolution(int& nump, double& freqstep);

    void SetT( double dtemperature);

    void SetZPLGamma(double ig);
    void UpdateFlagAcc(int iflag);
    
    void ReadCfun(string inp_file);


private:
    
private:
    
	double temperature;  
	double epsilonomega;

	int flagaccessories;
	// flagaccessories == 0 :
	// a3 function (correlation function) is calculated by FFT
	// a2 and a1 function is obtained by direct integration (includes omega = 0)
	// gfunction is calculated from accessories
	// flagaccessories == 1 :
	// all accessory functions are calculated by FFT
	// gfunction is calculated from accessories
	// flagaccessories == 2 :
	// a3 function (correlation function) is calculated by FFT
	// a2 and a1 function is obtained by direct integration (includes omega = 0)
	// gfunction is calculated by numerical integration of a3
  

private:
	asymptoticLF<double>* spectral_density;
	asymptoticLF<double> M_BAR_im;// only dataset for imaginary part in Fourier meaning
	asymptoticLF_complexv* lineshape;
	asymptoticLF_complexv* accessory_a1P;
	asymptoticLF_complexv* accessory_a2P;
	asymptoticLF_complexv* accessory_a3P;

    
	void Allocate_spectral_density();
	void Delete_spectral_density();

	void AllocateMBARFunction();
	void DeleteMBARFunction();
	void GenerateMBARFunction();
    
	void AllocateLineshape();
	void DeleteLineshape();
        
	void AllocateAccessoryFunction();
	void DeleteAccessoryFunction();

        void ConstructAccessoryFunction();
	void GenerateAccessoryArrays(
		double* dummySD, double* dummyPP, 
		complexv* funa1, complexv* funa2, complexv* funa3,
		int numt,double deltav);

	void CreateFDCorrelationSD(double* sfSSD,double* gfuSD,double deltav,int numt);
	void CreateSymmetricFDCorrelation(double* sfSSD,double* gfuSD,double deltav,int numt);
	void CreateLineshape();

};


