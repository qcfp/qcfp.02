

#ifndef CALCULATOR_1st_wf
#define	CALCULATOR_1st_wf


#pragma once

#include"../storage/storage.hpp"
#include"../complexv/complexv.hpp"
#include"../constructor-f-exciton/constructor-f-exciton.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../calculator_abs_secular_cumulant/calculator_abs_secular_cumulant.hpp"
#include"../averages_orient/averages_orient.hpp"
#include"../constants/constants.hpp"
#include"../interaction/interaction.hpp"
#include"../toolsFFT/toolsFFT.hpp"
#include"../propagatorExcitonWF/propagatorExcitonWF.hpp"

template<typename Tpropagator,typename Twavefunction>
class calculator_abs_wf:
public calculator_abs_secular_cumulant
{
    // generic function 
    // calculator based on the density matrix propagators
    // prepared as different projects
    // and used as hard-coded plugins
    
public:
    
    calculator_abs_wf<Tpropagator,Twavefunction>()
    :calculator_abs_secular_cumulant()
    {
        prop = 0;
        eigensys = false;
        outputstring = "#  calculator_abs_wf class signal\n";
        numEns = 1;
        propagationstep = 0.1;
        
    }
    calculator_abs_wf<Tpropagator,Twavefunction>(string& ifname)
    :calculator_abs_secular_cumulant()
    {
        prop = 0;
        eigensys = false;
        outputstring = "#  calculator_abs_wf class signal\n";
        numEns = 1;
        propagationstep = 0.1;
        

        toolsIO tio;
        ifstream istr(ifname.c_str());
        if(!istr.is_open())
        {
            cout<<"Error: input file cannot be read\nquitting.\n";
            return;
        }
        
        
        
        // all necessary information is read from a single file
        
        // reading the whole input file
        tio.ReadWholeFile(ifname,inputfile,1);
        
        
        
        tio.StreamSkipTrailers(istr);
        string tstr; getline(istr,tstr);
        ReadSystem(istr);
        ReadBath(istr);
        ReadExperimentLin(istr);
        
        // additional parameters for propagator
        string str = "";
        storage<double> arr(2);
        storage<int> iri(2);
        
        
        // timestep of propagation
        if(tio.LookUpAndReadSquareMatrix<double>(
            "MethodPropagatorNumericTimestep:",
            "Reading timestep of numerical propagation\n", 
            1,1, arr,inputfile.str()))
        {
            propagationstep = arr.data2D[0][0];
        }
        
        // size of thermal ensemble
        if(tio.LookUpAndReadSquareMatrix<int>(
            "MethodSizeOfThermalEnsemble:",
            "Reading numer of members in thermal ensemble\n", 
            1,1, iri,inputfile.str()))
        {
            numEns = iri.data2D[0][0];
        }
        
        
        
    }
    ~calculator_abs_wf<Tpropagator,Twavefunction>()
    {
        if(prop != 0)
            delete prop;
    }
    
    
    void Launch();
    
    double propagationstep;
    int numEns;
    
    
private:
    
    
    Twavefunction wfiLeft;
    Twavefunction wfiRight;
    
    Tpropagator* prop;
    
    
        
};



template<typename Tpropagator,typename Twavefunction>
void calculator_abs_wf<Tpropagator,Twavefunction>::
Launch()
{
    
    
//    if(!edips.IsAlloc())
//    {
//        MakeESSystem(1);
//        // this is necessary because the present calculations are done in eigenstate basis
//    }

	cout<<"# calculator_abs_wf::Launch initiated\n";

        if(prop == 0)
            prop = new Tpropagator(propagationstep);

        
    constants cst;
    double freSt = (ffre1-ifre1)/nump;
    double freF = ffre1-ifre1;
    // making time variables
    int timeN = nump/2;
    double timeS = cst.pi2/freF;
    double timeF = timeS*timeN;
    
    // system characteristics

    int& time1N = timeN;
    
    double time1i = 0;
    double& time1s = timeS;
            
    storage<complexv> dat3d(1);

    Twavefunction swf(*this);
    Twavefunction wfL0(*this);
    Twavefunction wfR0(*this);
    Twavefunction wfL1(*this);
    Twavefunction wfR1(*this);
    
    wfL0.meanEe = 0.5*(ffre1+ifre1);
    wfR0.meanEe = 0.5*(ffre1+ifre1);
    wfL1.meanEe = 0.5*(ffre1+ifre1);
    wfR1.meanEe = 0.5*(ffre1+ifre1);
    
    wfL0.meanEf = (ffre1+ifre1);
    wfR0.meanEf = (ffre1+ifre1);
    wfL1.meanEf = (ffre1+ifre1);
    wfR1.meanEf = (ffre1+ifre1);

    int numg = wfL0.numG;
    int nume = wfL0.numE;
    
    // serial code


	// now preparing transition properties (since calculation takes place in the site basis)

	storage<dvector3d> localedips(2);
	storage<dvector3d> localmdips(2);
	storage<dtensor3x3> localqdips(2);
	storage<int> skipper(2);
	int** skipflag=0;

	if(true)
	{
	    int numlev = 2;
	    if(bosonic) numlev = 3;
	    constructor_f_exciton exc_system(numlev,ham,dip);

		if(pos.IsAlloc())
			exc_system.AddDipolePositions(pos);
		if(d_m.IsAlloc())
			exc_system.AddMagneticDipoles(d_m);

		storage<dvector3d> ldip1(1);
		storage<dtensor3x3> lten1(1);
		storage<dvector3d> lmag1(1);

		exc_system.GetDips(ldip1);
		if(pos.IsAlloc())
		{
			exc_system.GetTens(lten1);
		}
		if(d_m.IsAlloc())
		{
			exc_system.GetMags(lmag1);
		}


		// making 2D matrix 
		localedips.Allocate(numg+nume,numg+nume);
		if(lten1.IsAlloc())
			localqdips.Allocate(numg+nume,numg+nume);
		if(lmag1.IsAlloc())
			localmdips.Allocate(numg+nume,numg+nume);


		    for(int ia = 0; ia<numE; ia++)
		    {
		        localedips.data2D[0][ia+1]=ldip1.data1D[ia];
		        localedips.data2D[ia+1][0]=ldip1.data1D[ia];

			if(lten1.IsAlloc())
			{
			        localqdips.data2D[0][ia+1]=lten1.data1D[ia];
			        localqdips.data2D[ia+1][0]=lten1.data1D[ia];
			}
			if(lmag1.IsAlloc())
			{
			        localmdips.data2D[0][ia+1]=lmag1.data1D[ia];
			        localmdips.data2D[ia+1][0]=lmag1.data1D[ia];
			}
		    }


		// transitions to skip:
		skipper.Allocate(numg+nume,numg+nume);
		skipflag=skipper.data2D;
		for(int ind=0;ind<numg+nume;ind++)
		for(int ins=0;ins<numg+nume;ins++)
		{
			if(localedips.data2D[ind][ins]==dvector3d(0.,0.,0.))
				skipflag[ind][ins] = 1;
		}
	}

        
    dat3d.Allocate(time1N*2);
        
        // running serial
    cout<<"Hello! I'm serial process\n";
      
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
	transition.PopulateSys(localedips.data2D,localqdips.data2D,localmdips.data2D);
        
    
    
    // looping over thermal ensemble
    for(int iens=0;iens<numEns; iens++)
    {
        
        // setting initial wavefunctions

        // next using only frequency grid
        swf.MakeThermal();
         
    // propagation of two branches and then taking the trace
    
    
    if(true)
    //     ig2 |  | ig2
    //       1-------
    //     ie4 |  | ig2
    //     ie1 |  | ig1
    //       0-------
    //     ig1 |  | ig1
    {
    	ios[0] = xomega1;
    	ios[1] = -xomega2;
	transition.PopulateFields(ies,iks,ios);

        wfL0.CopyState(swf);
        wfR0.CopyState(swf);

	// right line does nothing
        wfR1.CopyState(wfR0);

	// left line transition
        for(int ig1 = 0; ig1< numg; ig1++){
        for(int ie1 = 0; ie1< nume; ie1++){
	if(skipflag[ig1][ie1+numg]) continue;
	wfL1.CopyState(wfL0);
	wfL1.OpticalTransition01(ig1,ie1);
        transition.AssignLevels(0,ig1,ie1+numg); 

	for(int it1=0; it1<time1N; it1++){
        double time1 =  time1i + time1s*it1;
        prop->PropagateWF(time1,wfL1);
        prop->PropagateWF(time1,wfR1);

	complexv voverlap = wfL1.GetBathOverlap(wfR1);

	// left line transition
        for(int ig2 = 0; ig2< numg; ig2++){
        for(int ie4 = 0; ie4< nume; ie4++){
	if(skipflag[ie4+numg][ig2]) continue;
	wfL1.OpticalTransition10(ie4,ig2);
        transition.AssignLevels(1,ie4+numg,ig2); 

        complexv damplitude = transition.GetAveragedAmplitude();

        dat3d.data1D[it1] += coni*
            damplitude*wfL1.GetQuantumOverlap(wfR1)*voverlap*exp(-naturallinewidth*time1);
	// restoring
	wfL1.activemanifold = 1;
        }}
	}}}}

    // this is the end of general function
        
    }//looping over ensemble


	// now making FFT
	dat3d.data1D[0] /= 2.0;

    interpolationF<complexv> asymresFFTi(dat3d,0,timeS);
    interpolationF<complexv> asymresFFTf(0,(ffre1-ifre1)/nump,nump);
    
    toolsFFT fft;
    asymresFFTf = fft.executeN(asymresFFTi);
    fft.SwapSides(asymresFFTf);
    asymresFFTf.UpdateAxis(ifre1,(ffre1-ifre1)/nump);
    signal1d = asymresFFTf;


}






#endif	/* CALCULATOR_3RD_LEVELS_DMGEN_HPP */

