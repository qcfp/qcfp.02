

#ifndef CALCULATOR_3RD_wf
#define	CALCULATOR_3RD_wf


#pragma once

#include"../storage/storage.hpp"
#include"../complexv/complexv.hpp"
#include"../constructor-f-exciton/constructor-f-exciton.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../calculator_3rd_secular_cumulant/calculator_3rd_secular_cumulant.hpp"
#include"../averages_orient/averages_orient.hpp"
#include"../constants/constants.hpp"
#include"../propagatorExcitonWF/propagatorExcitonWF.hpp"

template<typename Tpropagator,typename Twavefunction>
class calculator_3rd_wf:
public calculator_3rd_secular_cumulant
{
    // generic third order response function 
    // calculator based on the density matrix propagators
    // prepared as different projects
    // and used as hard-coded plugins
    
    public:
        
        calculator_3rd_wf<Tpropagator,Twavefunction>():
            calculator_3rd_secular_cumulant()
        {
            prop = 0;
            eigensys = false;
            outputstring = "#  calculator_3rd_wf class signal\n";
            numEns = 1;
            propagationstep = 0.1;

        }       
        calculator_3rd_wf<Tpropagator,Twavefunction>(string& ifname):
            calculator_3rd_secular_cumulant()
        {
            prop = 0;
            eigensys = false;
            outputstring = "#  calculator_3rd_wf class signal\n";
            numEns = 1;
            propagationstep = 0.1;
            
            toolsIO tio;

            ifstream istr(ifname.c_str());
            if(!istr.is_open())
            {
                cout<<"Error: input file cannot be read\nquitting.\n";
                return;
            }
            // reading the whole input file
            tio.ReadWholeFile(ifname,inputfile,1);
            tio.StreamSkipTrailers(istr);
            string tstr; getline(istr,tstr);
            ReadSystem(istr);
            ReadBath(istr);
            ReadExperiment3rd(istr);

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
        
        ~calculator_3rd_wf<Tpropagator,Twavefunction>()
        {
            if(prop != 0)
                delete prop;
        }

        void LaunchGeneric(storage<complexv>& locdat3d, storage<complexv>& locdat3n, double* timeDP, int* timeIP);

	double propagationstep;
        int numEns;

        
    private:
        
        
        Twavefunction wfiLeft;
        Twavefunction wfiRight;
        
        Tpropagator* prop;
        
        
};



template<typename Tpropagator,typename Twavefunction>
void calculator_3rd_wf<Tpropagator,Twavefunction>::
LaunchGeneric(storage<complexv>& locdat3d, storage<complexv>& locdat3n, double* timeDP, int* timeIP)
{
    
    
	cout<<"################################################\n";
	cout<<"# calculator_3rd_wf::LaunchGeneric initiated\n";

        if(prop == 0)
            prop = new Tpropagator(propagationstep);

        
    // system characteristics

    int time3N = timeIP[2];
    int time2N = timeIP[1];
    int time1N = timeIP[0];
    
    double time3i = timeDP[5];
    double time2i = timeDP[4];
    double time1i = timeDP[3];
    double time3s = (timeDP[2]-timeDP[5])/time3N;
    double time2s = (timeDP[1]-timeDP[4])/time2N;
    double time1s = (timeDP[0]-timeDP[3])/time1N;
            
    storage<complexv> dat3d(3);
    storage<complexv> dat3n(3);

    Twavefunction swf(*this);
    Twavefunction wfL0(*this);
    Twavefunction wfR0(*this);
    Twavefunction wfL1(*this);
    Twavefunction wfR1(*this);
    Twavefunction wfL2(*this);
    Twavefunction wfR2(*this);
    Twavefunction wfL3(*this);
    Twavefunction wfR3(*this);
    
    int numg = wfL0.numG;

    // number of possible states in first manifold
    int nume = wfL0.numE;
    
    // number of possible states in second manifold
    int numf = wfL0.numF;

    int numeg= nume+numg;
    // serial code


	// now preparing transition properties (since calculation takes plave in the site basis)

	storage<dvector3d> localedips(2);
	storage<dvector3d> localmdips(2);
	storage<dtensor3x3> localqdips(2);
	storage<int> skipper(2);
	int** skipflag=0;

    
    // writing all transition properties in extended site basis
    // shifting the hamiltonian by the reprganization energies have been performedin the constructor
    // of wavefunctions
	if(true)
	{
	    int numlev = 2;
	    if(bosonic) numlev = 3;
	    constructor_f_exciton exc_system(numlev,ham,dip);
	    if(bosonic && dip2Corrections.IsAlloc()) 
	        exc_system.AddADipoles(dip2Corrections);

		if(pos.IsAlloc())
			exc_system.AddDipolePositions(pos);
		if(d_m.IsAlloc())
			exc_system.AddMagneticDipoles(d_m);

		storage<dvector3d> ldip1(1);
		storage<dvector3d> ldip2(2);
		storage<dtensor3x3> lten1(1);
		storage<dtensor3x3> lten2(2);
		storage<dvector3d> lmag1(1);
		storage<dvector3d> lmag2(2);

		exc_system.GetDips(ldip1);
		exc_system.GetDips2(ldip2);
		if(pos.IsAlloc())
		{
			exc_system.GetTens(lten1);
			exc_system.GetTens2(lten2);
		}
		if(d_m.IsAlloc())
		{
			exc_system.GetMags(lmag1);
			exc_system.GetMags2(lmag2);
		}


		// making 2D matrix 
		localedips.Allocate(numg+nume+numf,numg+nume+numf);
		if(lten1.IsAlloc())
			localqdips.Allocate(numg+nume+numf,numg+nume+numf);
		if(lmag1.IsAlloc())
			localmdips.Allocate(numg+nume+numf,numg+nume+numf);


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

			for(int ia = 0; ia<numF; ia++)
			{
			for(int ib = 0; ib<numE; ib++)
			{
			localedips.data2D[ia+1+numE][ib+1]=ldip2.data2D[ib][ia];
			localedips.data2D[ib+1][ia+1+numE]=ldip2.data2D[ib][ia];

			if(pos.IsAlloc())
			{
		            localqdips.data2D[ia+1+numE][ib+1]=lten2.data2D[ib][ia];
		            localqdips.data2D[ib+1][ia+1+numE]=lten2.data2D[ib][ia];
			}
			if(d_m.IsAlloc())
			{
		            localmdips.data2D[ia+1+numE][ib+1]=lmag2.data2D[ib][ia];
		            localmdips.data2D[ib+1][ia+1+numE]=lmag2.data2D[ib][ia];
			}

        		}
    			}

		// transitions to skip:
		skipper.Allocate(numg+nume+numf,numg+nume+numf);
		skipflag=skipper.data2D;
		for(int ind=0;ind<numg+nume+numf;ind++)
		for(int ins=0;ins<numg+nume+numf;ins++)
		{
			if(localedips.data2D[ind][ins]==dvector3d(0.,0.,0.))
				skipflag[ind][ins] = 1;
		}
	}

        
    dat3d.Allocate(timeIP[2],timeIP[1],timeIP[0]);
    dat3n.Allocate(timeIP[2],timeIP[1],timeIP[0]);
        
        // running serial
    cout<<"Hello! I'm serial process\n";
      
	// lab configuration
	dvector3d ies[4];
	dvector3d iks[4];
	double ios[4];
    
	// polarizations
	ies[0] = vecE1;
	ies[1] = vecE2;
	ies[2] = vecE3;
	ies[3] = vecE4;
    
	// wavevectors : these are laser directions (not of interactions) : arbitrary units - they are normalized later anyway
	iks[0] = veck1;
	iks[1] = veck2;
	iks[2] = veck3;
	iks[3] = veck4;

	// thse will be fourier frequencies . those of interactions +k = +w; -k = -w.
	ios[0] = xomega1;
	ios[1] = xomega2;
	ios[2] = xomega3;
	ios[3] = xomega4;

	int beyonddipole =0;
	if(eten.IsAlloc() || ed_m.IsAlloc()) beyonddipole =1;
	interaction transition(4, beyonddipole,averagingtype);
	transition.PopulateSys(localedips.data2D,localqdips.data2D,localmdips.data2D);
        
    
    
    // looping over thermal ensemble
    for(int iens=0;iens<numEns; iens++)
    {
        
        // setting initial wavefunctions

        // next using only frequency grid
        swf.MakeThermal();
         
    // propagation of two branches and then taking the trace
    
    
        /////////////////
    // K1 diagrams
    //if(false)
    
    if(coherentGSBK1 || transportGSBK1)
    // 1
    // coherent GSB K1
    //     ig6 |  | ig6
    //       3-------
    //     ie3 |  | ig6
    //     ie2 |  | ig5
    //       2-------
    //     ig3 |  | ig5
    //     ig2 |  | ig4
    //        -------1
    //     ig2 |  | ie4
    //     ig1 |  | ie1
    //        -------0
    //     ig1 |  | ig1
    
        
    {
    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

	std::cout<<"Running GSB\n";
        wfL0.CopyState(swf);
        wfR0.CopyState(swf);


	// left line nothing
        wfL1.CopyState(wfL0);

	// right line
        for(int ig1 = 0; ig1< numg; ig1++){
        for(int ie1 = 0; ie1< nume; ie1++){
	if(skipflag[ig1][ie1+numg]) continue;
        wfR1.CopyState(wfR0);
	wfR1.OpticalTransition01(ig1,ie1);
	transition.AssignLevels(0,ig1,ie1+numg); 

	for(int it1=0; it1<time1N; it1++){
        double time1 =  time1i + time1s*it1;

        prop->PropagateWF(time1,wfL1);
        prop->PropagateWF(time1,wfR1);

	// left line nothing
        wfL2.CopyState(wfL1);

	// right line
        for(int ig4 = 0; ig4< numg; ig4++){
        for(int ie4 = 0; ie4< nume; ie4++){
	if(skipflag[ie4+numg][ig4]) continue;
        wfR2.CopyState(wfR1);
	wfR2.OpticalTransition10(ie4,ig4);
        transition.AssignLevels(1,ie4+numg,ig4); 


	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;

        prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

	// right line nothing
        wfR3.CopyState(wfR2);

	// left line
        for(int ig3 = 0; ig3< numg; ig3++){
        for(int ie2 = 0; ie2< nume; ie2++){
	if(skipflag[ig3][ie2+numg]) continue;
        wfL3.CopyState(wfL2);
	wfL3.OpticalTransition01(ig3,ie2);
        transition.AssignLevels(2,ig3,ie2+numg); 


	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;

        prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // detection
        for(int ig6 = 0; ig6< numg; ig6++){
        for(int ie3 = 0; ie3< nume; ie3++){
	if(skipflag[ie3+numg][ig6]) continue;
	wfL3.OpticalTransition10(ie3,ig6);
        transition.AssignLevels(3,ie3+numg,ig6); 

        complexv damplitude = transition.GetAveragedAmplitude();
                        
        // finding electronic amplitude
        dat3d.data3D[it3][it2][it1] += cnni*
            damplitude*voverlap*wfL3.GetQuantumOverlap(wfR3);
	// restoration
	wfL3.activemanifold = 1;
        }}
	}}}}}}}}}}
      

       if(coherentESEK1 ||transportESEK1)
    // coherent ESE k1
    //     ig4 |  | ig4
    //       3-------
    //     ie6 |  | ig4
    //     ie5 |  | ig3
    //        -------2
    //     ie5 |  | ie3
    //     ie4 |  | ie2
    //       1-------
    //     ig2 |  | ie2
    //     ig1 |  | ie1
    //        -------0
    //     ig1 |  | ig1
    {
    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

	std::cout<<"Running ESE\n";
        wfL0.CopyState(swf);
        wfR0.CopyState(swf);

        // left line is nothing
	wfL1.CopyState(wfL0);

        // right line
        for(int ig1 = 0; ig1< numg; ig1++){
        for(int ie1 = 0; ie1< nume; ie1++){
	if(skipflag[ig1][ie1+numg]) continue;
	wfR1.CopyState(wfR0);
	wfR1.OpticalTransition01(ig1,ie1);
        transition.AssignLevels(0,ig1,ie1+numg); 

	for(int it1=0; it1<time1N; it1++){
        double time1 =  time1i + time1s*it1;
        prop->PropagateWF(time1,wfL1);
        prop->PropagateWF(time1,wfR1);

        // right line is nothing
	wfR2.CopyState(wfR1);

	// left line
        for(int ig2 = 0; ig2< numg; ig2++){
        for(int ie4 = 0; ie4< nume; ie4++){
	if(skipflag[ig2][ie4+numg]) continue;
	wfL2.CopyState(wfL1);
	wfL2.OpticalTransition01(ig2,ie4);
        transition.AssignLevels(1,ig2,ie4+numg); 
        
	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;
	prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

        // left line is nothing
	wfL3.CopyState(wfL2);

	// right line
        for(int ie3 = 0; ie3< nume; ie3++){
        for(int ig3 = 0; ig3< numg; ig3++){
	if(skipflag[ie3+numg][ig3]) continue;
	wfR3.CopyState(wfR2);
	wfR3.OpticalTransition10(ie3,ig3);
	transition.AssignLevels(2,ie3+numg,ig3); 

	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;
	prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // finding electronic amplitude
        for(int ig4 = 0; ig4< numg; ig4++){
        for(int ie6 = 0; ie6< nume; ie6++){
	if(skipflag[ie6+numg][ig4]) continue;
	wfL3.OpticalTransition10(ie6,ig4);
        transition.AssignLevels(3,ie6+numg,ig4); 

        complexv damplitude = transition.GetAveragedAmplitude();

        dat3d.data3D[it3][it2][it1] += cnni*
            damplitude*wfL2.GetQuantumOverlap(wfR3)*voverlap;

	wfL3.activemanifold = 1;
        }}
	}}}}}}}}}}

     

     
    if (coherentESAK1 || transportESAK1)
    // coherent ESA k1
    //     ie6 |  | ie6
    //       3-------
    //     if2 |  | ie6
    //     if1 |  | ie3
    //       2-------
    //     ie5 |  | ie3
    //     ie4 |  | ie2
    //       1-------
    //     ig2 |  | ie2
    //     ig1 |  | ie1
    //        -------0
    //     ig1 |  | ig1
    {
    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

        // left line
	std::cout<<"Running ESA\n";

        wfL0.CopyState(swf);
        wfR0.CopyState(swf);

        // left line does nothing
	wfL1.CopyState(wfL0);

	// right line transition
        for(int ig1 = 0; ig1< numg; ig1++){
        for(int ie1 = 0; ie1< nume; ie1++){
	if(skipflag[ig1][ie1+numg]) continue;
	wfR1.CopyState(wfR0);
	wfR1.OpticalTransition01(ig1,ie1);
        transition.AssignLevels(0,ig1,ie1+numg);  

	for(int it1=0; it1<time1N; it1++){
        double time1 =  time1i + time1s*it1;
        prop->PropagateWF(time1,wfL1);
        prop->PropagateWF(time1,wfR1);

	// right line does nothing
	wfR2.CopyState(wfR1);

	// left line transition
        for(int ig2 = 0; ig2< numg; ig2++){
        for(int ie4 = 0; ie4< nume; ie4++){
	if(skipflag[ig2][ie4+numg]) continue;
	wfL2.CopyState(wfL1);
	wfL2.OpticalTransition01(ig2,ie4);
        transition.AssignLevels(1,ig2,ie4+numg); 

	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;
	prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

	// right line does nothing
	wfR3.CopyState(wfR2);

	// left line transition to 2 exciton
        for(int ie5 = 0; ie5< nume; ie5++){
        for(int if1 = 0; if1< numf; if1++){
	if(skipflag[ie5+numg][if1+numeg]) continue;
	wfL3.CopyState(wfL2);
	wfL3.OpticalTransition12(ie5,if1);
        transition.AssignLevels(2,ie5+numg,if1+numeg); 

	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;
	prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // finding electronic amplitude
        for(int if2 = 0; if2< numf; if2++){
        for(int ie6 = 0; ie6< nume; ie6++){
	if(skipflag[if2+numeg][ie6+numg]) continue;
	wfL3.OpticalTransition21(if2,ie6);
        transition.AssignLevels(3,if2+numeg,ie6+numg); 

        complexv damplitude = transition.GetAveragedAmplitude();

        dat3n.data3D[it3][it2][it1] += cnni*
            damplitude*wfL3.GetQuantumOverlap(wfR1)*voverlap;
	// restoring
	wfL3.activemanifold = 2;
        }}
	}}}}}}}}}}




    /////////////////
    // K2 diagrams
    
    if(coherentGSBK2 || transportGSBK2)
     // coherent GSB K2
    //     ig6 |  | ig6
    //       3-------
    //     ie3 |  | ig6
    //     ie2 |  | ig5
    //       2-------
    //     ig3 |  | ig5
    //     ig4 |  | ig2
    //       1-------
    //     ie4 |  | ig2
    //     ie1 |  | ig1
    //       0-------
    //     ig1 |  | ig1
    {
    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

	std::cout<<"Running GSB\n";

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

	// right line does nothing
        wfR2.CopyState(wfR1);
        
	// left line transition
        for(int ig4 = 0; ig4< numg; ig4++){
        for(int ie4 = 0; ie4< nume; ie4++){
	if(skipflag[ie4+numg][ig4]) continue;
	wfL2.CopyState(wfL1);
	wfL2.OpticalTransition10(ie4,ig4);
        transition.AssignLevels(1,ie4+numg,ig4); 

	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;
	prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

	// right line does nothing
        wfR3.CopyState(wfR2);

	// left line transition
        for(int ig3 = 0; ig3< numg; ig3++){
        for(int ie2 = 0; ie2< nume; ie2++){
	if(skipflag[ig3][ie2+numg]) continue;
	wfL3.CopyState(wfL2);
	wfL3.OpticalTransition01(ig3,ie2);
        transition.AssignLevels(2,ig3,ie2+numg); 

	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;
	prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // finding electronic amplitude
        for(int ig6 = 0; ig6< numg; ig6++){
        for(int ie3 = 0; ie3< nume; ie3++){
	if(skipflag[ie3+numg][ig6]) continue;
	wfL3.OpticalTransition10(ie3,ig6); 
        transition.AssignLevels(3,ie3+numg,ig6); 
        // finding electronic amplitude
        complexv damplitude = transition.GetAveragedAmplitude();

        dat3d.data3D[it3][it2][it1] += cnni*
            damplitude*wfL3.GetQuantumOverlap(wfR3)*voverlap;
	// restoring
	wfL3.activemanifold = 1;
        }}
	}}}}}}}}}}

    
    if(coherentESEK2 || transportESEK2)
    // coherent ESE k2
    //     ig4 |  | ig4
    //        3------
    //     ie6 |  | ig4
    //     ie5 |  | ig3
    //        -------2
    //     ie5 |  | ie3
    //     ie2 |  | ie4
    //        -------1
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //       0-------
    //     ig1 |  | ig1
    {
    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

	std::cout<<"Running ESE\n";

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

	// left line does nothing
        wfL2.CopyState(wfL1);

	// right line transition
        for(int ig2 = 0; ig2< numg; ig2++){
        for(int ie4 = 0; ie4< nume; ie4++){
	if(skipflag[ig2][ie4+numg]) continue;	
	wfR2.CopyState(wfR1);
	wfR2.OpticalTransition01(ig2,ie4);
        transition.AssignLevels(1,ig2,ie4+numg); 

	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;
	prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

	// left line does nothing
        wfL3.CopyState(wfL2);

	// right line transition
        for(int ig3 = 0; ig3< numg; ig3++){
        for(int ie3 = 0; ie3< nume; ie3++){
	if(skipflag[ie3+numg][ig3]) continue;	
	wfR3.CopyState(wfR2);
	wfR3.OpticalTransition10(ie3,ig3); 
        transition.AssignLevels(2,ie3+numg,ig3); 

	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;
	prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // finding electronic amplitude
        for(int ig4 = 0; ig4< numg; ig4++){
        for(int ie6 = 0; ie6< nume; ie6++){
	if(skipflag[ie6+numg][ig4]) continue;
	wfL3.OpticalTransition10(ie6,ig4); 
        transition.AssignLevels(3,ie6+numg,ig4); 

        complexv damplitude = transition.GetAveragedAmplitude();

        dat3d.data3D[it3][it2][it1] += cnni*
            damplitude*wfL3.GetQuantumOverlap(wfR3)*voverlap;
	// restoring
	wfL3.activemanifold = 1;
        }}
	}}}}}}}}}}


      
    if(coherentESAK2 || transportESAK2)
    // coherent ESA k2
    //     ie6 |  | ie6
    //       3-------
    //     if2 |  | ie6
    //     if1 |  | ie3
    //       2-------
    //     ie5 |  | ie3
    //     ie4 |  | ie2
    //        -------1
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //       0-------
    //     ig1 |  | ig1
    {
    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

	std::cout<<"Running ESA\n";
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

	// left line does nothing
        wfL2.CopyState(wfL1);

	// right line transition
        for(int ig2 = 0; ig2< numg; ig2++){
        for(int ie2 = 0; ie2< nume; ie2++){
	if(skipflag[ig2][ie2+numg]) continue;
	wfR2.CopyState(wfR0);
	wfR2.OpticalTransition01(ig2,ie2);
        transition.AssignLevels(1,ig2,ie2+numg);

	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;
	prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

	// right line does nothing
        wfR3.CopyState(wfR2);

	// left line transition
        for(int ie5 = 0; ie5< nume; ie5++){
        for(int if1 = 0; if1< numf; if1++){
	if(skipflag[ie5+numg][if1+numeg]) continue;
	wfL3.OpticalTransition12(ie5,if1); 
	wfL3.CopyState(wfL2);
        transition.AssignLevels(2,ie5+numg,if1+numeg);  

	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;
	prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // finding electronic amplitude
        for(int ie6 = 0; ie6< nume; ie6++){
        for(int if2 = 0; if2< numf; if2++){
	if(skipflag[if2+numeg][ie6+numg]) continue;
	wfL3.OpticalTransition21(if2,ie6);
        transition.AssignLevels(3,if2+numeg,ie6+numg); 

        complexv damplitude = transition.GetAveragedAmplitude();

        dat3n.data3D[it3][it2][it1] += cnni*
            damplitude*wfL3.GetQuantumOverlap(wfR3)*voverlap;
	// restoring
	wfL3.activemanifold = 2;
        }}
	}}}}}}}}}}

    
    //////// now KIII
    
    if(coherentES1K3)
    // coherent ESA1 k3
    //     ig4 |  | ig4
    //        3-------
    //     ie4 |  | ig4
    //     ie3 |  | ig3
    //        2-------
    //     if2 |  | ig3
    //     if1 |  | ig2
    //        1-------
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //        0-------
    //     ig1 |  | ig1
    {
    	ios[0] = xomega1;
    	ios[1] = xomega2;
    	ios[2] = -xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

	std::cout<<"Running ES1\n";
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

	// right line does nothing
        wfR2.CopyState(wfR1);

	// left line transition
        for(int ie2 = 0; ie2< nume; ie2++){
        for(int if1 = 0; if1< numf; if1++){
	if(skipflag[ie2+numg][if1+numeg]) continue;
	wfL2.CopyState(wfL1);
	wfL2.OpticalTransition12(ie2,if1); 
        transition.AssignLevels(1,ie2+numg,if1+numeg); 

	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;
	prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

	// right line does nothing
        wfR3.CopyState(wfR2);

	// left line transition        
        for(int ie3 = 0; ie3< nume; ie3++){
        for(int if2 = 0; if2< numf; if2++){
	if(skipflag[if2+numeg][ie3+numg]) continue;
	wfL3.CopyState(wfL2);
	wfL3.OpticalTransition21(if2,ie3); 
        transition.AssignLevels(2,if2+numeg,ie3+numg); 

	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;
	prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // finding electronic amplitude
        for(int ig4 = 0; ig4< numg; ig4++){
        for(int ie4 = 0; ie4< nume; ie4++){
	if(skipflag[ie4+numg][ig4]) continue;
	wfL3.OpticalTransition10(ie4,ig4); 
        transition.AssignLevels(3,ie4+numg,ig4); 

        complexv damplitude = transition.GetAveragedAmplitude();

        dat3d.data3D[it3][it2][it1] += cnni*
            damplitude*wfL3.GetQuantumOverlap(wfR3)*voverlap;
	// restoring
	wfL3.activemanifold = 1;
        }}
	}}}}}}}}}}
    

    if(coherentES2K3)
    // coherent ESA2 k3
    //     ie4 |  | ie4
    //        3-------
    //     if3 |  | ie4
    //     if2 |  | ie3
    //        -------2
    //     if2 |  | ig3
    //     if1 |  | ig2
    //        1-------
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //        0-------
    //     ig1 |  | ig1
    {
    	ios[0] = xomega1;
    	ios[1] = xomega2;
    	ios[2] = -xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

	std::cout<<"Running ES2\n";
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

	// right line does nothing
        wfR2.CopyState(wfR1);

	// left line transition
        for(int ie2 = 0; ie2< nume; ie2++){
        for(int if1 = 0; if1< numf; if1++){
	if(skipflag[ie2+numg][if1+numeg]) continue;
	wfL2.CopyState(wfL1);
	wfL2.OpticalTransition12(ie2,if1); 
        transition.AssignLevels(1,ie2+numg,if1+numeg); 

	for(int it2=0; it2<time2N; it2++){
        double time2 =  time2i + time2s*it2;
	prop->PropagateWF(time1+time2,wfL2);
        prop->PropagateWF(time1+time2,wfR2);

	// left line does nothing
        wfL3.CopyState(wfL2);

	// right line transition
        for(int ig3 = 0; ig3< numg; ig3++){
        for(int ie3 = 0; ie3< nume; ie3++){
	if(skipflag[ig3][ie3+numg]) continue;
	wfR3.CopyState(wfR2);
	wfR3.OpticalTransition01(ig3,ie3);
        transition.AssignLevels(2,ig3,ie3+numg);

 	for(int it3=0; it3<time3N; it3++){
        double time3 =  time3i + time3s*it3;
	prop->PropagateWF(time1+time2+time3,wfL3);
        prop->PropagateWF(time1+time2+time3,wfR3);

	complexv voverlap = wfL3.GetBathOverlap(wfR3);

        // last electronic transition
        for(int ie4 = 0; ie4< nume; ie4++){
        for(int if3 = 0; if3< numf; if3++){
	if(skipflag[if3+numeg][ie4+numg]) continue;
	wfL2.OpticalTransition21(if3,ie4); 
        transition.AssignLevels(3,if3+numeg,ie4+numg); 

        complexv damplitude = transition.GetAveragedAmplitude();
        dat3n.data3D[it3][it2][it1] += cnni*
            damplitude*wfL3.GetQuantumOverlap(wfR3)*voverlap;
	// restoring
	wfL3.activemanifold = 2;
        }}
	}}}}}}}}}}
    

   
    int num1,num2,num3;
    dat3d.GetSize(num1,num2,num3);
    for(int it3=0; it3<num3; it3++)
        for(int it2=0; it2<num2; it2++)
            for(int it1=0; it1<num1; it1++)
            {
                locdat3d.data3D[it3][it2][it1] = cnni*dat3d.data3D[it3][it2][it1];
                locdat3n.data3D[it3][it2][it1] = cnni*dat3n.data3D[it3][it2][it1];
            }
  
    // this is the end of general function
        
    }//looping over ensemble

}






#endif	/* CALCULATOR_3RD_LEVELS_DMGEN_HPP */

