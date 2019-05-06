#pragma once
#include"../dvector3d/dvector3d.hpp"
#include"../complexv/complexv.hpp"
#include"../interaction/interaction.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
#include"../communicator_3rd/communicator_3rd.hpp"

#include"../toolsFFT/toolsFFT.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../storage/storage.hpp"
#include"../constants/constants.hpp"
#include"../numericalSD/numericalSD.hpp"
#include"../calculator_redfield/calculator_redfield.hpp"
#include"../propagatorExciton/propagatorExciton.hpp"

#include"../communicator_3rd/communicator_3rd.hpp"
#include"../calculator_abs_secular_cumulant/calculator_abs_secular_cumulant.hpp"



// no allocation / deallocation occurs here
// all pointers are passed from outside
template<typename Type>
class calculator_abs_reduced :
public calculator_abs_secular_cumulant
{
public:

    calculator_abs_reduced<Type>();
    calculator_abs_reduced<Type>(string& ifname);
    ~calculator_abs_reduced<Type>();

    // this function uses all input as pointers not allocating anything

    void Launch();


//    void ReadSpecific(ifstream& ifs);
	// nonsecular is defined in communicator
	int lindblad;
	double internaltimestep; // not good since this should go into propagatorExciton
  int internaltimepoints; // not good since this should go into propagatorExciton
	int flagMemoryWithCfun;

    private:

	// five propagators for  blocks

        Type* prop10;
};

template<typename Type>
calculator_abs_reduced<Type>::calculator_abs_reduced()
:calculator_abs_secular_cumulant()
{

	nonsecular = 0;
	markovian = 1;
	lindblad = 0;
  eigensys = false;
	flagMemoryWithCfun = 0;
	internaltimestep = 0.0;
	internaltimepoints = 0;

    prop10 = 0;

    outputstring = "#  calculator_abs_reduced class signal\n";

}
template<typename Type>
calculator_abs_reduced<Type>::calculator_abs_reduced(string& ifname)
:calculator_abs_secular_cumulant()
{
	storage<double> arr(2);
	storage<int> iri(2);

	nonsecular = 0;
	markovian = 1;
	lindblad = 0;
    eigensys = false;

    prop10 = 0;

    outputstring = "#  calculator_abs_reduced class signal\n";



    toolsIO tio;
    tio.ReadWholeFile(ifname,inputfile,1);

    // reading files
    ifstream istr(ifname.c_str());
    if(!istr.is_open())
    {
        cout<<"Error: input file cannot be read\nquitting.\n";
        return;
    }

    // reading the whole input file
    tio.StreamSkipTrailers(istr);
    string tstr; getline(istr,tstr);
    ReadSystem(istr);
    ReadBath(istr);
    ReadExperimentLin(istr);


	// additional parameters for propagator

    std::size_t found;
    found = experiment_complexity.find("Secular");
    if(found!=std::string::npos)
	nonsecular = 0;
    found = experiment_complexity.find("NonSecular");
    if(found!=std::string::npos)
	nonsecular = 1;

    found = experiment_complexity.find("Markovian");
    if(found!=std::string::npos)
	markovian=1;
    found = experiment_complexity.find("NonMarkovian");
    if(found!=std::string::npos)
		{
			markovian=0;
			flagMemoryWithCfun=0;
			nonsecular=1;
		}
		found = experiment_complexity.find("NonMarkovianCFun");
    if(found!=std::string::npos)
		{
			markovian=0;
			flagMemoryWithCfun=1;
			nonsecular=1;
		}
    found = experiment_complexity.find("Lindblad");
    if(found!=std::string::npos)
		{
			lindblad = 1;
			markovian=1;
			nonsecular = 1;
		}



		//implementation of internaltimeS and internaltimeN
	  internaltimestep=0;
		if(tio.LookUpAndReadSquareMatrix<double>(
	    			"MethodInternalTimeStep:",
	    			"Reading the timestep for correlation function\n",
	    			1,1, arr,inputfile.str()))
	  		{
	          internaltimestep = arr.data2D[0][0];
	     	 }

	  internaltimepoints=0;
		if(tio.LookUpAndReadSquareMatrix<int>(
	    			"MethodInternalTimePoints:",
	    			"Reading the number of points of correlation function\n",
	    			1,1, iri,inputfile.str()))
	  		{
	          internaltimepoints = iri.data2D[0][0];
	      	}


	////////////////////////////////////////////////////////////////////

}
template<typename Type>
calculator_abs_reduced<Type>::~calculator_abs_reduced()
{
    if(prop10 != 0)
        delete prop10;
}



template<typename Type>
void calculator_abs_reduced<Type>::Launch()
{
    //if(ready == 0)
    //{
    //    cout<<"Error: calculator_abs_secular_cumulant::Launch cannot proceed\n";
    //}
    cout<<"# Lauching calculator_abs_reduced::Launch\n";


    constants cst;
    double cfre = 0.5*(ifre1+ffre1);

    double freSt = (ffre1-ifre1)/nump;
    double freF = ffre1-ifre1;

    // making time variables
    int timeN = nump/2;
    double timeS = cst.pi2/freF;
    double timeF = timeS*timeN;

    storage<double> times1(1);
    if(timeN>1)
        times1.FillLinear(0.0,timeS,timeN);
    else if(timeN==1)
        times1.FillLinear(0.0,constants::smallepsilon,2);
    else
    {
        cout<<"Error: number of time points is incorrectly specified\n";
    }

    storage<complexv> dm1t10(3);
    storage<complexv> dat3d(1);
    dat3d.Allocate(nump);

    //making response functions (these will be used only in cumulant case
    //complexv* idat;
    //interpolationF<complexv> asymrest(0,timeS,timeN);
    //interpolationF<complexv> asymres(0,timeS,timeN);

    // system characteristics
    // number of possible states in ground manifold
    int& numg = numG;

    // number of possible states in first manifold
    int& nume = numE;

	// creating eigenstates (if missing)
	if(!evals.IsAlloc())
		MakeESSystem();

    FinalizeESSystem(0);

	// lab configuration
	dvector3d ies[2];
	dvector3d iks[2];
	double ios[2];

	// polarizations
	ies[0] = vecE1;
	ies[1] = vecE2;

	// wavevectors : these are laser directions (not of interactions) : arbitrary units - they are normalized later anyway
	iks[0] = veck1;
	iks[1] = veck2;

	// thse will be fourier frequencies . those of interactions +k = +w; -k = -w.
	ios[0] = xomega1;
	ios[1] = xomega2;

	int beyonddipole =0;
	if(eten.IsAlloc() || ed_m.IsAlloc()) beyonddipole =1;
	interaction transition(2, beyonddipole,averagingtype);
	transition.PopulateSys(edips.data2D,eten.data2D,ed_m.data2D);


	// Feinman diagrams
        if(true){
        if(prop10 == 0)
	{
          prop10 = new Type(*this);
	    		prop10->flagNonsecular = nonsecular;
					prop10->flagMarkovian=markovian;
					prop10->calcR->flagLindblad = lindblad;
					prop10->internaltimeS = internaltimestep;
    			prop10->internaltimeN =internaltimepoints;
					prop10->flagMemoryWithCfun = flagMemoryWithCfun;
					prop10->manifold0End = 0;

					prop10->SetBlock(10);
	}

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
			transition.PopulateFields(ies,iks,ios);


    //     ig2 |  | ig2
    //        -------
    //     ie4 |  | ig2
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1

        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
					//int ie1 = 0;

	    			transition.AssignLevels(0,ig1,ie1+numg);
            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

        for(int ig2 = 0; ig2< numg; ig2++)
        for(int ie4 = 0; ie4< nume; ie4++)
        {
	    			transition.AssignLevels(1,ie4+numg,ig2);

	    			complexv damplitude = transition.GetAveragedAmplitude();
            // time loops
            for(int it1=0; it1<timeN; it1++)
            {
            dat3d.data1D[it1] += coni*damplitude
                    *grPops.data1D[ig1]
//                    *exp( (coni*(cfre-12000.0)-20.0)*times1.data1D[it1]);
										*dm1t10.data3D[it1][ie4][ig2]*exp( (coni*cfre-naturallinewidth)*times1.data1D[it1]);
            }
        }}}



	// the rest is for time domain

    //cout<<"Error2\n";

    // shifting zero point due to Fourier transformation
    dat3d.data1D[0] /= 2.0;

    // doing Fourier
    interpolationF<complexv> asymresFFTi(dat3d,0,timeS);
    interpolationF<complexv> asymresFFTf(0,(ffre1-ifre1)/nump,nump);

    toolsFFT fft;
    //asymresFFTf = asymresFFTi;
    asymresFFTf = fft.executeN(asymresFFTi);

    //for(int ind = 0; ind< nump; ind++)
    //{
    //    cout<<"sig: "<<(ffre+ifre)/2+(ffre-ifre)/nump*ind<<"\t"<<asymresFFTf.Get((ffre-ifre)/nump*ind)<<"\n";
    //}



    //    cout<<asymresFFTf.Get(0)<<"\n";
    //    cout<<asymresFFTf.Get(0.1)<<"\n";

    // shifting frequencies
    fft.SwapSides(asymresFFTf);
    asymresFFTf.UpdateAxis(ifre1,(ffre1-ifre1)/nump);

    //for(int ind = 0; ind< nump; ind++)
    //{
    //    cout<<"sig: "<<ifre+(ffre-ifre)/nump*ind<<"\t"<<asymresFFTf.Get(ifre+(ffre-ifre)/nump*ind)<<"\n";
    //}

    //    cout<<asymresFFTf.Get(0)<<"\n";
    //    cout<<asymresFFTf.Get(0.1)<<"\n";

    // output
    signal1d = asymresFFTf;
}
