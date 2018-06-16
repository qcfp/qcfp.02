#pragma once


#include<string>
#include<iostream>
#include<fstream>
#include"../toolsFFT/toolsFFT.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../storage/storage.hpp"
#include"../constants/constants.hpp"
#include"../feinman2sideddiagram1/feinman2sideddiagram1.hpp"
#include"../interpolationF2d/interpolationF2d.hpp"
#include"../feinman2sideddiagram3/feinman2sideddiagram3.hpp"
#include"../propagatorM/propagatorM.hpp"
#include"../numericalSD/numericalSD.hpp"
#include"../calculator_redfield/calculator_redfield.hpp"
#include"../interaction/interaction.hpp"
#include"../dtensor3x3/dtensor3x3.hpp"
#include"../dvector3d/dvector3d.hpp"

using std::string;
using std::ifstream;
using std::ofstream;


#include"../communicator_3rd/communicator_3rd.hpp"

class calculator_3rd_secular_cumulant :
public communicator_3rd
{
    
public:
    storage<complexv> storedat1;
    storage<complexv> storedat2;
    storage<complexv> storedat3;

    double speedupsmallness;
    int noComplexLifetimes;
    int nonmarkovian;
public:

    calculator_3rd_secular_cumulant();
    calculator_3rd_secular_cumulant(string& ifname);

//    calculator_3rd_secular_cumulant(
//        storage<double>& ielevs,
//        storage<dvector3d>& iedips,
//        storage<complexv>& idephasings);

//    calculator_3rd_secular_cumulant(
//        storage<double>& ielevs,
//        storage<dvector3d>& iedips,
//        storage<complexv>& idephasings,
//        storage<double>& igrPops);
    
//    void ReadFile(string&);
//    void o_ReadFile(string& ifilename);

//   void SetupManifolds(int inumG, int inumE, int inumF);
    
    //void AddPhaseFunctions(storage<asymptoticLF_complexv>& imfun,storage<double>& iampls);
    //void AddLineshapes(storage<asymptoticLF_complexv>& igfun,storage<double>& iampls);
//    void AddLineshapes(storage< asymptoticLF_complexv>& igfun,
//                       storage< asymptoticLF_complexv>& imfun,
//                       storage<double>& igij,storage<double>& ikij,
//                       double itempr);

    //void AddDephasings();

    
    
    virtual interpolationF2d<complexv> LaunchTWW(double iifre2, double iffre2, double iifre3, double iffre3, int inump, double itf1);
    virtual interpolationF2d<complexv> LaunchWTW(double iifre1, double iffre1, double iifre3, double iffre3, int inump, double itf2);
    virtual interpolationF2d<complexv> LaunchWWT(double iifre1, double iffre1, double iifre2, double iffre2, int inump, double itf3);
    

    virtual void LaunchTTT();

    virtual interpolationF2d<complexv> LaunchTWW();
    // depends on ifre2, ffre2, ifre3, ffre3, nump, tf1

    virtual interpolationF2d<complexv> LaunchWTW();
    // depends on ifre1, ffre1, ifre3, ffre3, nump, tf2
    
    virtual interpolationF2d<complexv> LaunchWWT();
    // depends on ifre1, ffre1, ifre2, ffre2, nump, tf3
    
    virtual void LaunchBasic();
    
    virtual void PublishSpecial4(ofstream& ofs);
    
//    virtual void ReadSpecific(ifstream& ifs);
    

    // to launch generic time domain simulation
    virtual void  LaunchGeneric(storage<complexv>& locdat3d, storage<complexv>& locdat3n, double* timeDP, int* timeIP)
{
    


	cout<<"################################################\n";
	cout<<"# calculator_3rd_secular_cumulant::LaunchGeneric initiated\n";

    if(!edips.IsAlloc())
        MakeESSystem();
        // this is necessary because the present calculations are done in eigenstate basis

    FinalizeESSystem(nonmarkovian);

// 
    
    // system characteristics
    // number of possible states in ground manifold
    int& numg = numG;
    
    // number of possible states in first manifold
    int& nume = numE;
    
    // number of possible states in second manifold
    int& numf = numF;
    
    // local dipoles
    //dvector3d ids[4];

    
    int states[14];
    complexv freqs[3];
    
    // propagations only when transport is necessary
    // transport is only during t2
    //int tn2 = 0;
    //if( timeDP[1] >0 )
    //    tn2 = timeDP[1]/transportstep;
    //tn2+=2;
    //storage<double> trace2GG(3);
    //storage<double> trace2EE(3);
    
    propagatorM propgg;
    propagatorM propee;
	storage<double> propinigg(1);
	storage<double> propiniee(1);

         	propinigg.Allocate(numG);
        	propgg.InitMaster(ratesg);
         	propiniee.Allocate(numE);
        	propee.InitMaster(ratese);



    cout<<"# Initial population greens function matrix: GG:\n";
	for(int ig2 = 0; ig2< numG; ig2++)
       for(int ig1 = 0; ig1< numG; ig1++)
    {
    		propinigg.data1D[ig1]=1.0;
    		propgg.prepareGetAt(propinigg,ig2);
    		propinigg.data1D[ig1]=0.0;
    		double temp = 0.0;
    		cout<<propgg.GetAt(temp);
            if(ig1 == numG-1) cout<<"\n";
            else  cout<<"\t";
    }
    cout<<"# Final population greens function matrix: GG:\n";
	for(int ig2 = 0; ig2< numG; ig2++)
       for(int ig1 = 0; ig1< numG; ig1++)
    {
    		propinigg.data1D[ig1]=1.0;
    		propgg.prepareGetAt(propinigg,ig2);
    		propinigg.data1D[ig1]=0.0;
    		double temp = 0.0;
    		cout<<propgg.GetAt(timeDP[1]);
            if(ig1 == numG-1) cout<<"\n";
            else  cout<<"\t";
    }

    cout<<"# Initial population greens function matrix: EE:\n";
	for(int ig2 = 0; ig2< numE; ig2++)
       for(int ig1 = 0; ig1< numE; ig1++)
    {
    		propiniee.data1D[ig1]=1.0;
    		propee.prepareGetAt(propiniee,ig2);
    		propiniee.data1D[ig1]=0.0;
    		double temp = 0.0;
    		cout<<propee.GetAt(temp);
            if(ig1 == numE-1) cout<<"\n";
            else  cout<<"\t";
    }
    cout<<"# Final population greens function matrix: EE:\n";
	for(int ig2 = 0; ig2< numE; ig2++)
       for(int ig1 = 0; ig1< numE; ig1++)
    {
    		propiniee.data1D[ig1]=1.0;
    		propee.prepareGetAt(propiniee,ig2);
    		propiniee.data1D[ig1]=0.0;
    		double temp = 0.0;
    		cout<<propee.GetAt(timeDP[1]);
            if(ig1 == numE-1) cout<<"\n";
            else  cout<<"\t";
    }


    storage<complexv> dat3d(3);
    storage<complexv> dat3n(3);

    
    ///////////////////////////////////////////////////////////
    // MPI parallelization
#ifdef MPIPROCESSING
    
    dat3d.Allocate(1,1,1);
    dat3n.Allocate(1,1,1);
    
    int procid = 0;
    int numprocs = 0;
    /* find out MY process ID, and how many processes were started. */
    
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    
    
    // MPI - IIIIII  here the frequency loops are taken outside
    int time3N = timeIP[2];
    int time2N = timeIP[1];
    int time1N = timeIP[0];
    
    double time3i = timeDP[5];
    double time2i = timeDP[4];
    double time1i = timeDP[3];
    double time3s = (timeDP[2]-timeDP[5])/time3N;
    double time2s = (timeDP[1]-timeDP[4])/time2N;
    double time1s = (timeDP[0]-timeDP[3])/time1N;
    
    // setting up processes
    int loopLength =time1N*time2N*time3N;
    int loopSmall = loopLength / numprocs;
    if(loopLength % numprocs != 0)
        loopSmall += 1;
    int indini = procid*loopSmall;
    int indfin = (procid+1)*loopSmall;
    if(indfin>loopLength)
        indfin =loopLength;
    
    
    
    for ( int ii = 0; ii < numprocs; ++ii ) {
        if ( procid == ii ) {
            // my turn to write to the file
            // running parallel
            printf("Hello! I'm process %i out of %i processes\n", procid, numprocs);
            
            cout<<"my info: "<<procid<<" "<<loopLength<<" "<<indini<<" "<<indfin<<" "<<loopSmall<<"\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    
    double* mpidataexchageRd = new double[loopLength];
    double* mpidataexchageId = new double[loopLength];
    double* mpidataexchageRn = new double[loopLength];
    double* mpidataexchageIn = new double[loopLength];

    for(int ind = indini; ind<indfin; ind++)
    {
        int rem = ind;
        int itime1 = rem/ ( time2N*time3N );
        rem  = rem% ( time2N*time3N );
        int itime2 = rem/ ( time3N );
        rem  = rem% ( time3N );
        int itime3 = rem;
        
        double litT[6];
        int litN[3];
        
        litT[3] = litT[0] =  time1i+time1s * itime1;
        litT[4] = litT[1] =  time2i+time2s * itime2;
        litT[5] = litT[2] =  time3i+time3s * itime3;
        
        litN[0] = 1;
        litN[1] = 1;
        litN[2] = 1;
        
        dat3d.data3D[0][0][0] = 0.0;
        dat3n.data3D[0][0][0] = 0.0;
        
        
        //cout<<itime3<<" "<<itime2<<" "<<itime1<<"\n";
        
        ///// next are normal calculations of diagrams
        
#else
    // serial code
        
    dat3d.Allocate(timeIP[2],timeIP[1],timeIP[0]);
    dat3n.Allocate(timeIP[2],timeIP[1],timeIP[0]);
        
    double* litT = timeDP;
    int* litN = timeIP;
        
        
        // running serial
    cout<<"Hello! I'm serial process\n";

#endif

        
        // optimization
        int nume_numg =nume+numg;
        dvector3d** edipsdata2D = edips.data2D;
        double* evalsdata1D = evals.data1D;
        complexv** dephasingsdata2D=dephasings.data2D;
        asymptoticLF_complexv* gfundata1D = gfun.data1D;
        double*** gijdata3D = gij.data3D;
        int gfunGetSize = gfun.GetSize();
        double* grPopsdata1D = grPops.data1D;
        nstates = numG + numE + numF;
        int llll;
        
        
        if(noComplexLifetimes)
        {
            int n2,n1;
            dephasings.GetSize(n2,n1);
                for(int j1=0; j1<n2; j1++)
                    for(int j0=0; j0<n1;j0++)
                        dephasingsdata2D[j1][j0] = dephasingsdata2D[j1][j0].real();
                // notice that this overwrites dephasings.data2D as well
        }


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
	transition.PopulateSys(edips.data2D,eten.data2D,ed_m.data2D);
        
        /////////////////
    // K1 diagrams
    
    
    // 1
    // coherent GSB K1
    //     ig2 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie2 |  | ig2
    //        -------
    //     ig1 |  | ig2
    //     ig1 |  | ig2
    //        -------
    //     ig1 |  | ie1
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentGSBK1){

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

        for(int ig1 = 0; ig1< numg; ig1++)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ig2 = 0; ig2< numg; ig2++)
                    for(int ie2 = numg; ie2< nume_numg; ie2++)
                    {
                        
                        if(transport)
                            if(ig1 == ig2)
                                continue; // skipping

			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig2);
			transition.AssignLevels(2,ig1,ie2);
			transition.AssignLevels(3,ie2,ig2);
                        //ids[0] = edipsdata2D[ig1][ie1];
                        //ids[1] = edipsdata2D[ie1][ig2];
                        //ids[2] = edipsdata2D[ig1][ie2];
                        //ids[3] = edipsdata2D[ie2][ig2];
                        freqs[0] = (evalsdata1D[ig1]-evalsdata1D[ie1]-coni*dephasingsdata2D[ig1][ie1]);
                        freqs[1] = (evalsdata1D[ig1]-evalsdata1D[ig2]-coni*dephasingsdata2D[ig1][ig2]);
                        freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*dephasingsdata2D[ie2][ig2]);
                        
                        
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        
			if(nonmarkovian){
                        // making states indices
                        states[4 ] = ig2; /*|*/ states[4 ] = ig2; 
                        //--------------------------------
                        states[10] = ie2; /*|*/ states[11] = ig2; 
                        states[3 ] = ie2; /*|*/ states[5 ] = ig2; 
                        //--------------------------------
                        states[9 ] = ig1; /*|*/ states[12] = ig2; 
                        states[2 ] = ig1; /*|*/ states[6 ] = ig2; 
                        //--------------------------------
                        states[8 ] = ig1; /*|*/ states[13] = ie1; 
                        states[1 ] = ig1; /*|*/ states[7 ] = ie1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        

                        llll = 4;

                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
                        }
                        diagram.propagate(dat3d,litT,litN);
                        
                    }
	}

    // 1
    // transport GSB K1
    //     ig2 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie2 |  | ig2
    //        -------
    //     ig2 |  | ig2
    //     ig1 |  | ig1
    //        -------
    //     ig1 |  | ie1
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    
                    //cout<<"I am here\n";
    //if(false)
    if(transportGSBK1){

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

        if(transport)
            for(int ig1 = 0; ig1< numg; ig1++)
                if(grPopsdata1D[ig1]>speedupsmallness)
                for(int ig2 = 0; ig2< numg; ig2++)
                {
                    
                //    cout<<"I am here\n";
                    // making transport function
                    //asymptoticLF<double> propF(trace2GG.data3D[ig1][ig2],0,transportstep,tn2);
                    // now I have the full propagation function
                	propinigg.data1D[ig1]=1.0;
                	propgg.prepareGetAt(propinigg,ig2);
                	propinigg.data1D[ig1]=0.0;
                    
                    for(int ie1 = numg; ie1< nume_numg; ie1++)
                        for(int ie2 = numg; ie2< nume_numg; ie2++)
                        {
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,ig2,ie2);
			transition.AssignLevels(3,ie2,ig2);
                            //ids[0] = edipsdata2D[ig1][ie1];
                            //ids[1] = edipsdata2D[ie1][ig1];
                            //ids[2] = edipsdata2D[ig2][ie2];
                            //ids[3] = edipsdata2D[ie2][ig2];
                            freqs[0] = (evalsdata1D[ig1]-evalsdata1D[ie1]-coni*dephasingsdata2D[ig1][ie1]);
                            freqs[1] = 0.0;
                            freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*dephasingsdata2D[ie2][ig2]);
                            
                            
                            feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                            //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                            diagram.feinman2sideddiagram3smallparameter = speedupsmallness;
                            
			if(nonmarkovian){
                            // making states indices
                            
                            // transport GSB
                            states[4 ] = ig2; /*|*/ states[4 ] = ig2; 
                            //--------------------------------
                            states[10] = ie2; /*|*/ states[11] = ig2; 
                            states[3 ] = ie2; /*|*/ states[5 ] = ig2; 
                            //--------------------------------
                            states[9 ] = ig2; /*|*/ states[12] = ig2; 
                            states[2 ] = ig1; /*|*/ states[6 ] = ig1; 
                            //--------------------------------
                            states[8 ] = ig1; /*|*/ states[13] = ie1; 
                            states[1 ] = ig1; /*|*/ states[7 ] = ie1; 
                            //--------------------------------
                            states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                            
                            

                            llll = 10;
                            
                            diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
                            }
                            diagram.assignPropagator(&propgg);
                            
                            //diagram.multiplier(grPopsdata1D[ig1]);
                            
                            
                            //cout<<dat3d.data3D<<" "<<timeDP[0]<<" "<<timeDP[1]<<" "<<timeDP[2]<<" "<<timeIP[0]<<"\n";
                        diagram.propagate(dat3d,litT,litN);
                            
                        }
                }
    
    	}
    
    
    
    // 1
    // coherent ESE
    //     ig2 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie2 |  | ig2
    //        -------
    //     ie2 |  | ie1
    //     ie2 |  | ie1
    //        -------
    //     ig1 |  | ie1
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentESEK1){

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

        for(int ig1 = 0; ig1< numg; ig1++)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ig2 = 0; ig2< numg; ig2++)
                    for(int ie2 = numg; ie2< nume_numg; ie2++)
                    {
                        
                        if(transport)
                            if(ie1 == ie2)
                                continue; // skipping
                        
                        
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie2,ig1);
			transition.AssignLevels(2,ig2,ie1);
			transition.AssignLevels(3,ie2,ig2);
                        //ids[0] = edipsdata2D[ig1][ie1];
                        //ids[1] = edipsdata2D[ie2][ig1];
                        //ids[2] = edipsdata2D[ig2][ie1];
                        //ids[3] = edipsdata2D[ie2][ig2];
                        freqs[0] = (evalsdata1D[ig1]-evalsdata1D[ie1]-coni*dephasingsdata2D[ig1][ie1]);
                        freqs[1] = (evalsdata1D[ie2]-evalsdata1D[ie1]-coni*dephasingsdata2D[ie2][ie1]);
                        freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*dephasingsdata2D[ie2][ig2]);
                        
                        
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        
                        
			if(nonmarkovian){
                        // making states indices
                        
                        states[4 ] = ig2; /*|*/ states[4 ] = ig2; 
                        //--------------------------------
                        states[10] = ie2; /*|*/ states[11] = ig2; 
                        states[3 ] = ie2; /*|*/ states[5 ] = ig2; 
                        //--------------------------------
                        states[9 ] = ie2; /*|*/ states[12] = ie1; 
                        states[2 ] = ie2; /*|*/ states[6 ] = ie1; 
                        //--------------------------------
                        states[8 ] = ig1; /*|*/ states[13] = ie1; 
                        states[1 ] = ig1; /*|*/ states[7 ] = ie1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        
                        
                        //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);

                        llll = 6;
                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}                        
                        //diagram.multiplier(grPopsdata1D[ig1]);
                        
                        diagram.propagate(dat3d,litT,litN);
                    }
    
    	}
    
    
    // 1
    // transport ESE
    //     ig2 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie2 |  | ig2
    //        -------
    //     ie2 |  | ie2
    //     ie1 |  | ie1
    //        -------
    //     ig1 |  | ie1
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    
    if(transportESEK1){

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

        if(transport)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                {
                    
                    //asymptoticLF<double> propF(trace2EE.data3D[ie1-numg][ie2-numg],0,transportstep,tn2);
                    
                	propiniee.data1D[ie1-numg]=1.0;
                	propee.prepareGetAt(propiniee,ie2-numg);
                	propiniee.data1D[ie1-numg]=0.0;

                	for(int ig1 = 0; ig1< numg; ig1++)
                        if(grPopsdata1D[ig1]>speedupsmallness)
                        for(int ig2 = 0; ig2< numg; ig2++)
                        {
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,ig2,ie2);
			transition.AssignLevels(3,ie2,ig2);
//                            ids[0] = edipsdata2D[ig1][ie1];
//                            ids[1] = edipsdata2D[ie1][ig1];
//                            ids[2] = edipsdata2D[ig2][ie2];
//                            ids[3] = edipsdata2D[ie2][ig2];
                            freqs[0] = (evalsdata1D[ig1]-evalsdata1D[ie1]-coni*dephasingsdata2D[ig1][ie1]);
                            freqs[1] = 0.0;
                            freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*dephasingsdata2D[ie2][ig2]);
                            
                            
                            feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);                            
                            diagram.feinman2sideddiagram3smallparameter = speedupsmallness;

			if(nonmarkovian){
                            // making states indices
                            
                            // transport GSB
                            states[4 ] = ig2; /*|*/ states[4 ] = ig2; 
                            //--------------------------------
                            states[10] = ie2; /*|*/ states[11] = ig2; 
                            states[3 ] = ie2; /*|*/ states[5 ] = ig2; 
                            //--------------------------------
                            states[9 ] = ie2; /*|*/ states[12] = ie2; 
                            states[2 ] = ie1; /*|*/ states[6 ] = ie1; 
                            //--------------------------------
                            states[8 ] = ig1; /*|*/ states[13] = ie1; 
                            states[1 ] = ig1; /*|*/ states[7 ] = ie1; 
                            //--------------------------------
                            states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                            
                            
                            //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
  
                            llll = 10;
                            diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}
                            diagram.assignPropagator(&propee);
                            
                            //diagram.multiplier(grPopsdata1D[ig1]);
                            
                            diagram.propagate(dat3d,litT,litN);
                            
                        }
                }
    
    
    	}
    
    // 1
    // coherent ESA
    //     ie1 |  | ie1
    //        -------
    //     if1 |  | ie1
    //     if1 |  | ie1
    //        -------
    //     ie2 |  | ie1
    //     ie2 |  | ie1
    //        -------
    //     ig1 |  | ie1
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentESAK1){

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

        for(int ig1 = 0; ig1< numg; ig1++)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                    for(int if1 = nume_numg; if1< (numf+nume_numg); if1++)
                    {
                        
                        if(transport)
                            if(ie1 == ie2)
                                continue; // skipping
                        
                        
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie2,ig1);
			transition.AssignLevels(2,if1,ie2);
			transition.AssignLevels(3,ie1,if1);
//                        ids[0] = edipsdata2D[ig1][ie1];
//                        ids[1] = edipsdata2D[ie2][ig1];
//                        ids[2] = edipsdata2D[if1][ie2];
//                        ids[3] = edipsdata2D[ie1][if1];
                        freqs[0] = (evalsdata1D[ig1]-evalsdata1D[ie1]-coni*dephasingsdata2D[ig1][ie1]);
                        freqs[1] = (evalsdata1D[ie2]-evalsdata1D[ie1]-coni*dephasingsdata2D[ie2][ie1]);
                        freqs[2] = (evalsdata1D[if1]-evalsdata1D[ie1]-coni*dephasingsdata2D[if1][ie1]);
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        
                        // making states indices
                        
			if(nonmarkovian){
                        states[4 ] = ie1; /*|*/ states[4 ] = ie1; 
                        //--------------------------------
                        states[10] = if1; /*|*/ states[11] = ie1; 
                        states[3 ] = if1; /*|*/ states[5 ] = ie1; 
                        //--------------------------------
                        states[9 ] = ie2; /*|*/ states[12] = ie1; 
                        states[2 ] = ie2; /*|*/ states[6 ] = ie1; 
                        //--------------------------------
                        states[8 ] = ig1; /*|*/ states[13] = ie1; 
                        states[1 ] = ig1; /*|*/ states[7 ] = ie1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        
                        
                        
//                        cout<<"dephasings1 = "<<dephasingsdata2D[ig1][ie1]<<"\n";
//                        cout<<"dephasings2 = "<<dephasingsdata2D[ie2][ie1]<<"\n";
//                        cout<<"dephasings3 = "<<dephasingsdata2D[if1][ie1]<<"\n";
                        
                        //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);

                        llll = 2;
                        
                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}                        
                        //diagram.multiplier(grPopsdata1D[ig1]);
                        
                        diagram.propagate(dat3n,litT,litN);
                        
                    }
    
    
    	}
    
    
    // 1
    // transport ESA
    //     ie2 |  | ie2
    //        -------
    //     if1 |  | ie2
    //     if1 |  | ie2
    //        -------
    //     ie2 |  | ie2
    //     ie1 |  | ie1
    //        -------
    //     ig1 |  | ie1
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
//        if(false)
    if(transportESAK1){

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

        if(transport)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                {
                    //asymptoticLF<double> propF(trace2EE.data3D[ie1-numg][ie2-numg],0,transportstep,tn2);
                    
                	propiniee.data1D[ie1-numg]=1.0;
                	propee.prepareGetAt(propiniee,ie2-numg);
                	propiniee.data1D[ie1-numg]=0.0;
                    
                    for(int ig1 = 0; ig1< numg; ig1++)
                        if(grPopsdata1D[ig1]>speedupsmallness)
                        for(int if1 = nume_numg; if1< (nume_numg+numf); if1++)
                        {
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,if1,ie2);
			transition.AssignLevels(3,ie2,if1);
//                            ids[0] = edipsdata2D[ig1][ie1];
//                            ids[1] = edipsdata2D[ie1][ig1];
//                            ids[2] = edipsdata2D[if1][ie2];
//                            ids[3] = edipsdata2D[ie2][if1];
                            freqs[0] = (evalsdata1D[ig1]-evalsdata1D[ie1]-coni*dephasingsdata2D[ig1][ie1]);
                            freqs[1] = 0.0;
                            freqs[2] = (evalsdata1D[if1]-evalsdata1D[ie2]-coni*dephasingsdata2D[if1][ie2]);
                            
                            
                            feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                            //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                            diagram.feinman2sideddiagram3smallparameter = speedupsmallness;
                            
			if(nonmarkovian){
                            // making states indices
                            
                            // transport ESA
                            states[4 ] = ie2; /*|*/ states[4 ] = ie2; 
                            //--------------------------------
                            states[10] = if1; /*|*/ states[11] = ie2; 
                            states[3 ] = if1; /*|*/ states[5 ] = ie2; 
                            //--------------------------------
                            states[9 ] = ie2; /*|*/ states[12] = ie2; 
                            states[2 ] = ie1; /*|*/ states[6 ] = ie1; 
                            //--------------------------------
                            states[8 ] = ig1; /*|*/ states[13] = ie1; 
                            states[1 ] = ig1; /*|*/ states[7 ] = ie1; 
                            //--------------------------------
                            states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                            
                            
                            
                            llll = 10;
                            
                            diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}
                            diagram.assignPropagator(&propee);
                            
                            //diagram.multiplier(grPopsdata1D[ig1]);
                            
                            diagram.propagate(dat3n,litT,litN);
                            
                            //cout<<"inside trESAK1\n";
                            
                        }
                }
    
    	}
    // all K1 diagrams have been calculated
    
    /////////////////
    // K2 diagrams
    
    
    // 1
    // coherent GSB K2
    //     ig1 |  | ig1
    //        -------
    //     ie2 |  | ig1
    //     ie2 |  | ig1
    //        -------
    //     ig2 |  | ig1
    //     ig2 |  | ig1
    //        -------
    //     ie1 |  | ig1
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentGSBK2){

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


        for(int ig1 = 0; ig1< numg; ig1++)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ig2 = 0; ig2< numg; ig2++)
                    for(int ie2 = numg; ie2< nume_numg; ie2++)
                    {
                        
                        if(transport)
                            if(ig1 == ig2)
                                continue; // skipping

			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig2);
			transition.AssignLevels(2,ig2,ie2);
			transition.AssignLevels(3,ie2,ig1);
//                        ids[0] = edipsdata2D[ig1][ie1];
//                        ids[1] = edipsdata2D[ie1][ig2];
//                        ids[2] = edipsdata2D[ig2][ie2];
//                        ids[3] = edipsdata2D[ie2][ig1];
                        freqs[0] = (evalsdata1D[ie1]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie1][ig1]);
                        freqs[1] = (evalsdata1D[ig2]-evalsdata1D[ig1]-coni*dephasingsdata2D[ig2][ig1]);
                        freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie2][ig1]);
                        
                        //freqs[0] = -cfre1 + (evalsdata1D[ig1]-evalsdata1D[ie1]-coni*100.);
                        //freqs[1] = (evalsdata1D[ig1]-evalsdata1D[ig2]-coni*dephasingsdata2D[ig1][ig2]);
                        //freqs[2] = -cfre3 + (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*100.);
                        
                        //cout<<"freqs:"<<freqs[0]<<"\n";
                        
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        
			if(nonmarkovian){
                        // making states indices
                        states[4 ] = ig1; /*|*/ states[4 ] = ig1; 
                        //--------------------------------
                        states[10] = ie2; /*|*/ states[11] = ig1; 
                        states[3 ] = ie2; /*|*/ states[5 ] = ig1; 
                        //--------------------------------
                        states[9 ] = ig2; /*|*/ states[12] = ig1; 
                        states[2 ] = ig2; /*|*/ states[6 ] = ig1; 
                        //--------------------------------
                        states[8 ] = ie1; /*|*/ states[13] = ig1; 
                        states[1 ] = ie1; /*|*/ states[7 ] = ig1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        
                        
                        //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                        
                        llll = 1;
                        
                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}                        
                        //diagram.multiplier(grPopsdata1D[ig1]);
                        
                        diagram.propagate(dat3d,litT,litN);
                        
                    }
    
    	}
    //    for(int i2 = 0; i2<timeN; i2++)
    //    for(int i1 = 0; i1<timeN; i1++)
    //        cout<<"GSB:"<<dat3d.data3D[i2][0][i1]<<"\n";
    
    
    
    // 1
    // transport GSB K2
    //     ig2 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie2 |  | ig2
    //        -------
    //     ig2 |  | ig2
    //     ig1 |  | ig1
    //        ------
    //     ie1 |  | ig1
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    //if(false)
    if(transportGSBK2){

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


        if(transport)
            for(int ig1 = 0; ig1< numg; ig1++)
                if(grPopsdata1D[ig1]>speedupsmallness)
                for(int ig2 = 0; ig2< numg; ig2++)
                {
                    // making transport function
                    //asymptoticLF<double> propF(trace2GG.data3D[ig1][ig2],0,transportstep,tn2);
                    // now I have the full propagation function
                	propinigg.data1D[ig1]=1.0;
                	propgg.prepareGetAt(propinigg,ig2);
                	propinigg.data1D[ig1]=0.0;
                    
                    for(int ie1 = numg; ie1< nume_numg; ie1++)
                        for(int ie2 = numg; ie2< nume_numg; ie2++)
                        {
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,ig2,ie2);
			transition.AssignLevels(3,ie2,ig2);
//                            ids[0] = edipsdata2D[ig1][ie1];
//                            ids[1] = edipsdata2D[ie1][ig1];
//                            ids[2] = edipsdata2D[ig2][ie2];
//                            ids[3] = edipsdata2D[ie2][ig2];
                            freqs[0] = (evalsdata1D[ie1]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie1][ig1]);
                            freqs[1] = 0.0;
                            freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*dephasingsdata2D[ie2][ig2]);
                            
                            
                            feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                            
                            //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                            diagram.feinman2sideddiagram3smallparameter = speedupsmallness;
                            
			if(nonmarkovian){
                            // making states indices
                            
                            // transport GSB
                            states[4 ] = ig2; /*|*/ states[4 ] = ig2; 
                            //--------------------------------
                            states[10] = ie2; /*|*/ states[11] = ig2; 
                            states[3 ] = ie2; /*|*/ states[5 ] = ig2; 
                            //--------------------------------
                            states[9 ] = ig2; /*|*/ states[12] = ig2; 
                            states[2 ] = ig1; /*|*/ states[6 ] = ig1; 
                            //--------------------------------
                            states[8 ] = ie1; /*|*/ states[13] = ig1; 
                            states[1 ] = ie1; /*|*/ states[7 ] = ig1; 
                            //--------------------------------
                            states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                            
                            
                            llll = 9;
                            
                            diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
                            
                            diagram.assignPropagator(&propgg);
                            
                            //diagram.multiplier(grPopsdata1D[ig1]);
			}                            
                            diagram.propagate(dat3d,litT,litN);
                            
                        }
                }
    
    
    
    	}
    
    
    // 1
    // coherent ESE
    //     ig2 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie2 |  | ig2
    //        -------
    //     ie2 |  | ie1
    //     ie2 |  | ie1
    //        -------
    //     ie2 |  | ig1
    //     ie2 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentESEK2){

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


        for(int ig1 = 0; ig1< numg; ig1++)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ig2 = 0; ig2< numg; ig2++)
                    for(int ie2 = numg; ie2< nume_numg; ie2++)
                    {
                        
                        if(transport)
                            if(ie1 == ie2)
                                continue; // skipping
                        
			transition.AssignLevels(0,ig1,ie2);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,ig2,ie1);
			transition.AssignLevels(3,ie2,ig2);                        
//                        ids[0] = edipsdata2D[ig1][ie2];
//                        ids[1] = edipsdata2D[ie1][ig1];
//                        ids[2] = edipsdata2D[ig2][ie1];
//                        ids[3] = edipsdata2D[ie2][ig2];
                        freqs[0] = (evalsdata1D[ie2]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie2][ig1]);
                        freqs[1] = (evalsdata1D[ie2]-evalsdata1D[ie1]-coni*dephasingsdata2D[ie2][ie1]);
                        freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*dephasingsdata2D[ie2][ig2]);
                        
                        
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                        
                        // making states indices
                        
			if(nonmarkovian){

                        states[4 ] = ig2; /*|*/ states[4 ] = ig2; 
                        //--------------------------------
                        states[10] = ie2; /*|*/ states[11] = ig2; 
                        states[3 ] = ie2; /*|*/ states[5 ] = ig2; 
                        //--------------------------------
                        states[9 ] = ie2; /*|*/ states[12] = ie1; 
                        states[2 ] = ie2; /*|*/ states[6 ] = ie1; 
                        //--------------------------------
                        states[8 ] = ie2; /*|*/ states[13] = ig1; 
                        states[1 ] = ie2; /*|*/ states[7 ] = ig1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        
                        
                        
                        llll = 7;
                        
                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}                        
                        //diagram.multiplier(grPopsdata1D[ig1]);
                        
                        diagram.propagate(dat3d,litT,litN);
                    }
    
    	}
    
    
    // 1
    // transport ESE
    //     ig2 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie2 |  | ig2
    //        -------
    //     ie2 |  | ie2
    //     ie1 |  | ie1
    //        -------
    //     ie1 |  | ig1
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    if(transportESEK2){

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


        if(transport)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                {
                    
                    //asymptoticLF<double> propF(trace2EE.data3D[ie1-numg][ie2-numg],0,transportstep,tn2);
                    
                	propiniee.data1D[ie1-numg]=1.0;
                	propee.prepareGetAt(propiniee,ie2-numg);
                	propiniee.data1D[ie1-numg]=0.0;

                	for(int ig1 = 0; ig1< numg; ig1++)
                        if(grPopsdata1D[ig1]>speedupsmallness)
                        for(int ig2 = 0; ig2< numg; ig2++)
                        {
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,ig2,ie2);
			transition.AssignLevels(3,ie2,ig2);
//                            ids[0] = edipsdata2D[ig1][ie1];
//                            ids[1] = edipsdata2D[ie1][ig1];
//                            ids[2] = edipsdata2D[ig2][ie2];
//                            ids[3] = edipsdata2D[ie2][ig2];
                            freqs[0] = (evalsdata1D[ie1]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie1][ig1]);
                            freqs[1] = 0.0;
                            freqs[2] = (evalsdata1D[ie2]-evalsdata1D[ig2]-coni*dephasingsdata2D[ie2][ig2]);
                            
                            
                            feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                            //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                            diagram.feinman2sideddiagram3smallparameter = speedupsmallness;
                            
                            
                            
			if(nonmarkovian){
                            // making states indices
                            
                            // transport GSB
                            states[4 ] = ig2; /*|*/ states[4 ] = ig2; 
                            //--------------------------------
                            states[10] = ie2; /*|*/ states[11] = ig2; 
                            states[3 ] = ie2; /*|*/ states[5 ] = ig2; 
                            //--------------------------------
                            states[9 ] = ie2; /*|*/ states[12] = ie2; 
                            states[2 ] = ie1; /*|*/ states[6 ] = ie1; 
                            //--------------------------------
                            states[8 ] = ie1; /*|*/ states[13] = ig1; 
                            states[1 ] = ie1; /*|*/ states[7 ] = ig1; 
                            //--------------------------------
                            states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                            
                            
                            llll = 9;
                            
                            diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}
                            diagram.assignPropagator(&propee);
                            
                            //diagram.multiplier(grPopsdata1D[ig1]);
                            
                            diagram.propagate(dat3d,litT,litN);
                            
                        }
                }
    
    
    	}
    
    // 1
    // coherent ESA
    //     ie1 |  | ie1
    //        -------
    //     if1 |  | ie1
    //     if1 |  | ie1
    //        -------
    //     ie2 |  | ie1
    //     ie2 |  | ie1
    //        -------
    //     ie2 |  | ig1
    //     ie2 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentESAK2){

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


        for(int ig1 = 0; ig1< numg; ig1++)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                    for(int if1 = nume_numg; if1< (numf+nume_numg); if1++)
                    {
                        
                        if(transport)
                            if(ie1 == ie2)
                                continue; // skipping
                        
			transition.AssignLevels(0,ig1,ie2);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,if1,ie2);
			transition.AssignLevels(3,ie1,if1);
//                        ids[0] = edipsdata2D[ig1][ie2];
//                        ids[1] = edipsdata2D[ie1][ig1];
//                        ids[2] = edipsdata2D[if1][ie2];
//                        ids[3] = edipsdata2D[ie1][if1];
                        freqs[0] = (evalsdata1D[ie2]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie2][ig1]);
                        freqs[1] = (evalsdata1D[ie2]-evalsdata1D[ie1]-coni*dephasingsdata2D[ie2][ie1]);
                        freqs[2] = (evalsdata1D[if1]-evalsdata1D[ie1]-coni*dephasingsdata2D[if1][ie1]);
                        
                        
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                        
                        
			if(nonmarkovian){
                        // making states indices
                        
                        states[4 ] = ie1; /*|*/ states[4 ] = ie1; 
                        //--------------------------------
                        states[10] = if1; /*|*/ states[11] = ie1; 
                        states[3 ] = if1; /*|*/ states[5 ] = ie1; 
                        //--------------------------------
                        states[9 ] = ie2; /*|*/ states[12] = ie1; 
                        states[2 ] = ie2; /*|*/ states[6 ] = ie1; 
                        //--------------------------------
                        states[8 ] = ie2; /*|*/ states[13] = ig1; 
                        states[1 ] = ie2; /*|*/ states[7 ] = ig1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        
                        
                        llll = 3;
                        
                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}                        
                        //diagram.multiplier(grPopsdata1D[ig1]);
                        
                        diagram.propagate(dat3n,litT,litN);
                        
                    }
    
    
    	}
    
    
    // 1
    // transport ESA
    //     ie2 |  | ie2
    //        -------
    //     if1 |  | ie2
    //     if1 |  | ie2
    //        -------
    //     ie2 |  | ie2
    //     ie1 |  | ie1
    //        -------
    //     ie1 |  | ig1
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    if(transportESAK2){

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


        if(transport)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                {
                    //asymptoticLF<double> propF(trace2EE.data3D[ie1-numg][ie2-numg],0,transportstep,tn2);
                    
                	propiniee.data1D[ie1-numg]=1.0;
                	propee.prepareGetAt(propiniee,ie2-numg);
                	propiniee.data1D[ie1-numg]=0.0;
                    
                    for(int ig1 = 0; ig1< numg; ig1++)
                    if(grPopsdata1D[ig1]>speedupsmallness)
                        for(int if1 = nume_numg; if1< (nume_numg+numf); if1++)
                        {
			transition.AssignLevels(0,ig1,ie1);
			transition.AssignLevels(1,ie1,ig1);
			transition.AssignLevels(2,if1,ie2);
			transition.AssignLevels(3,ie2,if1);
//                            ids[0] = edipsdata2D[ig1][ie1];
//                            ids[1] = edipsdata2D[ie1][ig1];
//                            ids[2] = edipsdata2D[if1][ie2];
//                            ids[3] = edipsdata2D[ie2][if1];
                            freqs[0] = (evalsdata1D[ie1]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie1][ig1]);
                            freqs[1] = 0.0;
                            freqs[2] = (evalsdata1D[if1]-evalsdata1D[ie2]-coni*dephasingsdata2D[if1][ie2]);
                            
                            
                            feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                            //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                            diagram.feinman2sideddiagram3smallparameter = speedupsmallness;
                          
                            
			if(nonmarkovian){
                            // making states indices
                            
                            // transport ESA
                            states[4 ] = ie2; /*|*/ states[4 ] = ie2; 
                            //--------------------------------
                            states[10] = if1; /*|*/ states[11] = ie2; 
                            states[3 ] = if1; /*|*/ states[5 ] = ie2; 
                            //--------------------------------
                            states[9 ] = ie2; /*|*/ states[12] = ie2; 
                            states[2 ] = ie1; /*|*/ states[6 ] = ie1; 
                            //--------------------------------
                            states[8 ] = ie1; /*|*/ states[13] = ig1; 
                            states[1 ] = ie1; /*|*/ states[7 ] = ig1; 
                            //--------------------------------
                            states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                            
                            
                            llll = 9;
                            
                            diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}
                            diagram.assignPropagator(&propee);
                            
                            //diagram.multiplier(grPopsdata1D[ig1]);
                            
                            diagram.propagate(dat3n,litT,litN);
                            
                        }
                }
    
    
    
    	}
    
    
    
    
    //////// now KIII
    
    // 1
    // coherent ESA
    //     ig1 |  | ig1
    //        -------
    //     ie1 |  | ig1
    //     ie1 |  | ig1
    //        -------
    //     if1 |  | ig1
    //     if1 |  | ig1
    //        -------
    //     ie2 |  | ig1
    //     ie2 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentES1K3){

    	ios[0] = xomega1;
    	ios[1] = xomega2;
    	ios[2] = -xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


        for(int ig1 = 0; ig1< numg; ig1++)
                if(grPopsdata1D[ig1]>speedupsmallness)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                    for(int if1 = nume_numg; if1< (numf+nume_numg); if1++)
                    {                        
                        
			transition.AssignLevels(0,ig1,ie2);
			transition.AssignLevels(1,ie2,if1);
			transition.AssignLevels(2,if1,ie1);
			transition.AssignLevels(3,ie1,ig1);
//                        ids[0] = edipsdata2D[ig1][ie2];
//                        ids[1] = edipsdata2D[ie2][if1];
//                        ids[2] = edipsdata2D[if1][ie1];
//                        ids[3] = edipsdata2D[ie1][ig1];
                        
                        freqs[0] = (evalsdata1D[ie2]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie2][ig1]);
                        freqs[1] = (evalsdata1D[if1]-evalsdata1D[ig1]-coni*dephasingsdata2D[if1][ig1]);
                        freqs[2] = (evalsdata1D[ie1]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie1][ig1]);
                        
                        
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);

			if(nonmarkovian){
                        // making states indices
                        
                        states[4 ] = ig1; /*|*/ states[4 ] = ig1; 
                        //--------------------------------
                        states[10] = ie1; /*|*/ states[11] = ig1; 
                        states[3 ] = ie1; /*|*/ states[5 ] = ig1; 
                        //--------------------------------
                        states[9 ] = if1; /*|*/ states[12] = ig1; 
                        states[2 ] = if1; /*|*/ states[6 ] = ig1; 
                        //--------------------------------
                        states[8 ] = ie2; /*|*/ states[13] = ig1; 
                        states[1 ] = ie2; /*|*/ states[7 ] = ig1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        
                      
                        llll = 1;
                        
                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}                        
                        //diagram.multiplier(grPopsdata1D[ig1]);
                        
                        diagram.propagate(dat3d,litT,litN);
                        
                    }
    
    	}
    // 1
    // coherent ESA
    //     ie1 |  | ie1
    //        -------
    //     if1 |  | ie1
    //     if1 |  | ie1
    //        -------
    //     if1 |  | ig1
    //     if1 |  | ig1
    //        -------
    //     ie2 |  | ig1
    //     ie2 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
    if(coherentES2K3){

    	ios[0] = xomega1;
    	ios[1] = xomega2;
    	ios[2] = -xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);



        for(int ig1 = 0; ig1< numg; ig1++)
                if(grPopsdata1D[ig1]>speedupsmallness)
            for(int ie1 = numg; ie1< nume_numg; ie1++)
                for(int ie2 = numg; ie2< nume_numg; ie2++)
                    for(int if1 = nume_numg; if1< (numf+nume_numg); if1++)
                    {
                        
			transition.AssignLevels(0,ig1,ie2);
			transition.AssignLevels(1,ie2,if1);
			transition.AssignLevels(2,ig1,ie1);
			transition.AssignLevels(3,ie1,if1);
//                        ids[0] = edipsdata2D[ig1][ie2];
//                        ids[1] = edipsdata2D[ie2][if1];
//                        ids[2] = edipsdata2D[ig1][ie1];
//                        ids[3] = edipsdata2D[ie1][if1];
                        freqs[0] = (evalsdata1D[ie2]-evalsdata1D[ig1]-coni*dephasingsdata2D[ie2][ig1]);
                        freqs[1] = (evalsdata1D[if1]-evalsdata1D[ig1]-coni*dephasingsdata2D[if1][ig1]);
                        freqs[2] = (evalsdata1D[if1]-evalsdata1D[ie1]-coni*dephasingsdata2D[if1][ie1]);
                        
                        
                        feinman2sideddiagram3 diagram( freqs, transition.GetAveragedAmplitude()*grPopsdata1D[ig1]);
                        //feinman2sideddiagram3 diagram( freqs, ids, ies, averagingtype);
                        
			if(nonmarkovian){
                        // making states indices
                        
                        states[4 ] = ie1; /*|*/ states[4 ] = ie1; 
                        //--------------------------------
                        states[10] = if1; /*|*/ states[11] = ie1; 
                        states[3 ] = if1; /*|*/ states[5 ] = ie1; 
                        //--------------------------------
                        states[9 ] = if1; /*|*/ states[12] = ig1; 
                        states[2 ] = if1; /*|*/ states[6 ] = ig1; 
                        //--------------------------------
                        states[8 ] = ie2; /*|*/ states[13] = ig1; 
                        states[1 ] = ie2; /*|*/ states[7 ] = ig1; 
                        //--------------------------------
                        states[0 ] = ig1; /*|*/ states[0 ] = ig1; 
                        
                        

                        llll = 3;
                        
                        diagram.assignLineshape(llll,gfundata1D,gijdata3D,gfunGetSize,states, nstates);
			}                        
                        //diagram.multiplier(grPopsdata1D[ig1]);
                        
                        diagram.propagate(dat3n,litT,litN);
                        
                    }


	}
#ifdef MPIPROCESSING
// finalizing MPI processing
        
        // making the response function
        // assigning result
        dat3d.data3D[0][0][0] =  cnni*dat3d.data3D[0][0][0];
        dat3n.data3D[0][0][0] =  cnni*dat3n.data3D[0][0][0];
        mpidataexchageRd[ind-indini] = dat3d.data3D[0][0][0].real();
        mpidataexchageId[ind-indini] = dat3d.data3D[0][0][0].imag();
        mpidataexchageRn[ind-indini] = dat3n.data3D[0][0][0].real();
        mpidataexchageIn[ind-indini] = dat3n.data3D[0][0][0].imag();
        
        
    }
    
    
    
    
    for ( int ii = 0; ii < numprocs; ++ii ) {
        if ( procid == ii ) {
            // my turn to write to the file
            // running parallel
            printf("Process %i out of %i processes entering data exchange\n", procid, numprocs);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    
    // MPI - IIIIII data exchange
    if (procid != 0) {
        MPI_Send(&indini, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&indfin, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(mpidataexchageRd, indfin-indini, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        MPI_Send(mpidataexchageId, indfin-indini, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
        MPI_Send(mpidataexchageRn, indfin-indini, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
        MPI_Send(mpidataexchageIn, indfin-indini, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
        
        
        for(int it3=0; it3<time3N; it3++)
            for(int it2=0; it2<time2N; it2++)
                for(int it1=0; it1<time1N; it1++)
                {
                    locdat3d.data3D[it3][it2][it1] = 0.0;
                    locdat3n.data3D[it3][it2][it1] = 0.0;
                }

        
    } else {
        
        // do data collection  and result accumulation
        for(int pid = 0; pid<numprocs; pid++)
        {
            if(pid != 0)
            {
                
                MPI_Status status;
                
                MPI_Recv(&indini, 1, MPI_INT, pid, 1, MPI_COMM_WORLD,&status);
                MPI_Recv(&indfin, 1, MPI_INT, pid, 2, MPI_COMM_WORLD,&status);
                MPI_Recv(mpidataexchageRd, indfin-indini, MPI_DOUBLE, pid, 3, MPI_COMM_WORLD,&status);
                MPI_Recv(mpidataexchageId, indfin-indini, MPI_DOUBLE, pid, 4, MPI_COMM_WORLD,&status);
                MPI_Recv(mpidataexchageRn, indfin-indini, MPI_DOUBLE, pid, 5, MPI_COMM_WORLD,&status);
                MPI_Recv(mpidataexchageIn, indfin-indini, MPI_DOUBLE, pid, 6, MPI_COMM_WORLD,&status);
            }
            
            for(int ind = indini; ind<indfin; ind++)
            {
                int rem = ind;
                int itime1 = rem/ ( time2N*time3N );
                rem  = rem% ( time2N*time3N );
                int itime2 = rem/ ( time3N );
                rem  = rem% ( time3N );
                int itime3 = rem;
                
                locdat3d.data3D[itime3][itime2][itime1] = complexv(mpidataexchageRd[ind-indini],mpidataexchageId[ind-indini]);
                locdat3n.data3D[itime3][itime2][itime1] = complexv(mpidataexchageRn[ind-indini],mpidataexchageIn[ind-indini]);
            }
        }
    }
    delete[] mpidataexchageRd;
    delete[] mpidataexchageId;
    delete[] mpidataexchageRn;
    delete[] mpidataexchageIn;
    
    
    
    
    for ( int ii = 0; ii < numprocs; ++ii ) {
        if ( procid == ii ) {
            // my turn to write to the file
            // running parallel
            printf("Process %i out of %i processes finished assignment\n", procid, numprocs);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    
#else

    int num1,num2,num3;
    dat3d.GetSize(num1,num2,num3);
    for(int it3=0; it3<num3; it3++)
        for(int it2=0; it2<num2; it2++)
            for(int it1=0; it1<num1; it1++)
            {
                locdat3d.data3D[it3][it2][it1] = cnni*dat3d.data3D[it3][it2][it1];
                locdat3n.data3D[it3][it2][it1] = cnni*dat3n.data3D[it3][it2][it1];
            }
    
#endif
    
    // this is the end of general function
}




}; // end of class

