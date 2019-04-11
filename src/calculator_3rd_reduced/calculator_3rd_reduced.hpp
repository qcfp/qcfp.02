/* 
 * File:   calculator_3rd_reduced.hpp
 * Author: dariusa
 *
 * Created on November 5, 2015, 8:32 PM
 */
#pragma once

#include"../calculator_3rd_secular_cumulant/calculator_3rd_secular_cumulant.hpp"
#include"../averages_orient/averages_orient.hpp"
#include"../constants/constants.hpp"
#include"../interaction/interaction.hpp"
#include"../dtensor3x3/dtensor3x3.hpp"
#include"../dvector3d/dvector3d.hpp"

template<typename Type>
class calculator_3rd_reduced:
public calculator_3rd_secular_cumulant
{
    // generic third order response function 
    // calculator based on the density matrix propagators
    // prepared as different projects
    // and used as hard-coded plugins
    
    public:
        
        calculator_3rd_reduced<Type>():
            calculator_3rd_secular_cumulant()
        {
            nonsecular = 0;
            markovian = 1;
            lindblad = 0;

            prop00 = 0;
            prop10 = 0;
            prop11 = 0;
            prop20 = 0;
            prop21 = 0;
            eigensys = false;
            outputstring = "#  calculator_3rd_reduced class signal\n";
        }
        calculator_3rd_reduced<Type>(string& ifname);
        ~calculator_3rd_reduced<Type>();

        void LaunchGeneric(storage<complexv>& locdat3d, storage<complexv>& locdat3n, double* timeDP, int* timeIP);
    
        

	// nonsecular is defined in communicator
	int markovian;
	int lindblad;

        
    private:
        
	// five propagators for  blocks
        Type* prop00;
        Type* prop10;
        Type* prop11;
        Type* prop20;
        Type* prop21;
};




template<typename Type>
calculator_3rd_reduced<Type>::calculator_3rd_reduced(string& ifname):
    calculator_3rd_secular_cumulant()
{

        // now reading
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

    std::size_t found;
    
    nonsecular = 0;
    found = experiment_complexity.find("Secular");
    if(found!=std::string::npos)
	nonsecular = 0;
    found = experiment_complexity.find("NonSecular");
    if(found!=std::string::npos)
	nonsecular = 1;

    markovian = 1;
    found = experiment_complexity.find("Markovian");
    if(found!=std::string::npos)
	markovian=1;
    found = experiment_complexity.find("NonMarkovian");
    if(found!=std::string::npos)
	markovian=0;

    lindblad = 0;
    found = experiment_complexity.find("Lindblad");
    if(found!=std::string::npos)
	lindblad = 1;


    
}
template<typename Type>
calculator_3rd_reduced<Type>::~calculator_3rd_reduced()
{
    if(prop00 != 0)
        delete prop00;
    if(prop10 != 0)
        delete prop10;
    if(prop11 != 0)
        delete prop11;
    if(prop20 != 0)
        delete prop20;
    if(prop21 != 0)
        delete prop21;
}




template<typename Type>
void
calculator_3rd_reduced<Type>::LaunchGeneric
(storage<complexv>& locdat3d, storage<complexv>& locdat3n, double* timeDP, int* timeIP)
{
    

    cout<<"################################################\n";
    cout<<"# calculator_3rd_reduced::LaunchGeneric initiated\n";

        
    if(!edips.IsAlloc())
        MakeESSystem();
        // this is necessary because the present calculations are done in eigenstate basis
   
    FinalizeESSystem(nonmarkovian);


    // system characteristics
    // number of possible states in ground manifold
    int& numg = numG;
    
    // number of possible states in first manifold
    int& nume = numE;
    
    // number of possible states in second manifold
    int& numf = numF;
    
    // local dipoles
    //dvector3d ids[4];
    //dvector3d ies[4];
    
    
    //ies[0] = vecE1;
    //ies[1] = vecE2;
    //ies[2] = vecE3;
    //ies[3] = vecE4;
    
    //averages_orient obj_avor(averagingtype);


    int time3N = timeIP[2];
    int time2N = timeIP[1];
    int time1N = timeIP[0];
    
    double time3i = timeDP[5];
    double time2i = timeDP[4];
    double time1i = timeDP[3];
    double time3s = (timeDP[2]-timeDP[5])/time3N;
    double time2s = (timeDP[1]-timeDP[4])/time2N;
    double time1s = (timeDP[0]-timeDP[3])/time1N;
    
    storage<double> times1(1); 
    if(time1N>1)
        times1.FillLinear(time1i,time1s,time1N);
    else if(time1N==1)
        times1.FillLinear(time1i,constants::smallepsilon,2);
    else
    {
        cout<<"Error: number of time points is incorrectly specified\n";
    }
    storage<double> times2(1);    
    if(time2N>1)
        times2.FillLinear(time2i,time2s,time2N);
    else if(time2N==1)
        times2.FillLinear(time2i,constants::smallepsilon,2);
    else
    {
        cout<<"Error: number of time points is incorrectly specified\n";
    }
    storage<double> times3(1);    
    if(time3N>1)
        times3.FillLinear(time3i,time3s,time3N);
    else if(time3N==1)
        times3.FillLinear(time3i,constants::smallepsilon,2);
    else
    {
        cout<<"Error: number of time points is incorrectly specified\n";
    }
        
    storage<complexv> dat3d(3);
    storage<complexv> dat3n(3);

    
    // serial code
        
    dat3d.Allocate(timeIP[2],timeIP[1],timeIP[0]);
    dat3n.Allocate(timeIP[2],timeIP[1],timeIP[0]);
        
        // running serial
    cout<<"Hello! I'm serial process\n";
      
    // optimization
    int nume_numg =nume+numg;
    dvector3d** edipsdata2D = edips.data2D;
    double* grPopsdata1D = grPops.data1D;
        
        
    storage<complexv> dm1t00(3);
    storage<complexv> dm2t00(3);
    storage<complexv> dm3t00(3);

    storage<complexv> dm1t10(3);
    storage<complexv> dm2t10(3);
    storage<complexv> dm3t10(3);
    
    storage<complexv> dm1t20(3);
    storage<complexv> dm2t20(3);
    storage<complexv> dm3t20(3);
    
    storage<complexv> dm1t11(3);
    storage<complexv> dm2t11(3);
    storage<complexv> dm3t11(3);
    
    storage<complexv> dm1t21(3);
    storage<complexv> dm2t21(3);
    storage<complexv> dm3t21(3);

    // for non-markovian GetHistory
    //storage<complexv> dm2auxt10(2);

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
    
    if(coherentGSBK1 || transportGSBK1)
    { 
        if(prop00 == 0)
        {
            prop00 = new Type(*this);
	    	prop00->flagNonsecular = nonsecular;
            prop00->flagMarkovian=markovian;
            prop00->calcR->flagLindblad = lindblad;
            prop00->SetBlock(0);
        }

        if(prop10 == 0)
        {
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
            prop10->flagMarkovian=markovian;
            prop10->calcR->flagLindblad = lindblad;
            prop10->SetBlock(10);
        }

    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);


    
    // 1
    // coherent GSB K1
    //     ig6 |  | ig6
    //        -------
    //     ie3 |  | ig6
    //     ie2 |  | ig5
    //        -------
    //     ig3 |  | ig5
    //     ig2 |  | ig4
    //        -------
    //     ig2 |  | ie4
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
            transition.AssignLevels(0,ig1,ie1+numg);
            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

        for(int ig2 = 0; ig2< numg; ig2++)
        for(int ig4 = 0; ig4< numg; ig4++)
        {
            prop00->SetDeltaDM(ig2,ig4);
            dm2t00=prop00->PropagateDM(times2);
   
        for(int ig5 = 0; ig5< numg; ig5++)
        for(int ie2 = 0; ie2< nume; ie2++)
        {
            prop10->SetDeltaDM(ie2,ig5);
            dm3t10=prop10->PropagateDM(times3);

        for(int ie4 = 0; ie4< nume; ie4++)
        {
	    transition.AssignLevels(1,ie4+numg,ig4);
        for(int ig3 = 0; ig3< numg; ig3++)
        {
	    transition.AssignLevels(2,ig3,ie2+numg);
        for(int ie3 = 0; ie3< nume; ie3++)
        for(int ig6 = 0; ig6< numg; ig6++)
        {

	    transition.AssignLevels(3,ie3+numg,ig6);

        // result
//            double damplitude = obj_avor.rot_av_dip_4(ies,ids);
	    complexv damplitude = transition.GetAveragedAmplitude();
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3d.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *conj(dm1t10.data3D[it1][ie4][ig2])
                    *dm2t00.data3D[it2][ig3][ig5]
                    *dm3t10.data3D[it3][ie3][ig6];
                    
                    
 /*                   if(!markovian)
                    {
                        
    //     ig6 |  | ig6
    //        -------
    //     ie3 |  | ig6
    //     ie2 |  | ig5
    //        -------
    //     ig3 |  | ig5
    //     ig2 |  | ig4
    //        -------
    //     ig2 |  | ie4
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    

                        for(int iea =0; iea<numE; iea++)
                        for(int ieb =0; ieb<numE; ieb++)
                        for(int iec =0; iec<numE; iec++)
                        for(int ied =0; ied<numE; ied++)
                        {
                            

                        kernel10 = M(ie4,ie5, ie5,ie1) sum e5
                        
                        kernel20 = 0
                        kernel21 = 0;

                        kernel30 = M(ie3,ie5, ie5,ie2) sum e5
                        kernel31 = 0;
                        kernel32 = M(ie3,ie2, ie4,ie1)
                        
                        
                        prop10->SetDeltaDM(ie1,ig1);
                        
                    }
                    
                    
 */                   
                    
            }
        }}}}}}
    }

        
     
    
        if(coherentESEK1 ||transportESEK1)
    { 

        if(prop11 == 0)
	{
            prop11 = new Type(*this);
	    	prop11->flagNonsecular = nonsecular;
		prop11->flagMarkovian=markovian;
		prop11->calcR->flagLindblad = lindblad;
		prop11->SetBlock(11);
	}
        if(prop10 == 0)
	{
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
		prop10->flagMarkovian=markovian;
		prop10->calcR->flagLindblad = lindblad;
		prop10->SetBlock(10);
	}
    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

    // coherent ESE k1
    //     ig4 |  | ig4
    //        -------
    //     ie6 |  | ig4
    //     ie5 |  | ig3
    //        -------
    //     ie5 |  | ie3
    //     ie4 |  | ie2
    //        -------
    //     ig2 |  | ie2
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1
    
        for(int ie2 = 0; ie2< nume; ie2++)
        for(int ie4 = 0; ie4< nume; ie4++)
        {
            prop11->SetDeltaDM(ie4,ie2);
            dm2t11=prop11->PropagateDM(times2);
   
            //if(!markovian){
            //    prop10->PropagateDM(times2);
            //    dm2auxt10=prop10->GetHistory();
            //}
            
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
	    transition.AssignLevels(0,ie1+numg,ig1);
            //ids[0] = edipsdata2D[ig1][ie1+numg];

            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

        for(int ig3 = 0; ig3< numg; ig3++)
        for(int ie5 = 0; ie5< nume; ie5++)
        {
            prop10->SetDeltaDM(ie5,ig3);
            //if(!markovian){
            //    prop10->SetHistory(dm2auxt10);
            //}
            dm3t10=prop10->PropagateDM(times3);

        for(int ig2 = 0; ig2< numg; ig2++)
        {
	    transition.AssignLevels(1,ie4+numg,ig2);
            //ids[1] = edipsdata2D[ie4+numg][ig2];
        for(int ie3 = 0; ie3< nume; ie3++)
        {
	    transition.AssignLevels(2,ie3+numg,ig3);
            //ids[2] = edipsdata2D[ie3+numg][ig3];
        for(int ie6 = 0; ie6< nume; ie6++)
        for(int ig4 = 0; ig4< numg; ig4++)
        {

	    transition.AssignLevels(3,ie6+numg,ig4);
            //ids[3] = edipsdata2D[ie6+numg][ig4];

            // result
	    complexv damplitude = transition.GetAveragedAmplitude();
//            double damplitude = obj_avor.rot_av_dip_4(ies,ids);
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3d.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *conj(dm1t10.data3D[it1][ie2][ig2])
                    *dm2t11.data3D[it2][ie5][ie3]
                    *dm3t10.data3D[it3][ie6][ig4];
            }
        }}}}}}
    }
        
    
    
     
    if (coherentESAK1 || transportESAK1)
    { 

        if(prop10 == 0)
	{
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
		prop10->flagMarkovian=markovian;
		prop10->calcR->flagLindblad = lindblad;
		prop10->SetBlock(10);
	}
        if(prop11 == 0)
	{
            prop11 = new Type(*this);
	    	prop11->flagNonsecular = nonsecular;
		prop11->flagMarkovian=markovian;
		prop11->calcR->flagLindblad = lindblad;
		prop11->SetBlock(11);
	}
        if(prop21 == 0)
	{
            prop21 = new Type(*this);
	    	prop21->flagNonsecular = nonsecular;
		prop21->flagMarkovian=markovian;
		prop21->calcR->flagLindblad = lindblad;
		prop21->SetBlock(21);
	}
    	ios[0] = -xomega1;
    	ios[1] = xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

   // coherent ESA k1
    //     ie6 |  | ie6
    //        -------
    //     if2 |  | ie6
    //     if1 |  | ie3
    //        -------
    //     ie5 |  | ie3
    //     ie4 |  | ie2
    //        -------
    //     ig2 |  | ie2
    //     ig1 |  | ie1
    //        -------
    //     ig1 |  | ig1

	//int if1 = 0;
	//int ie3 = 0;
        for(int if1 = 0; if1< numf; if1++)
        for(int ie3 = 0; ie3< nume; ie3++)
        {
            prop21->SetDeltaDM(if1,ie3);
            dm3t21=prop21->PropagateDM(times3);

        for(int ie2 = 0; ie2< nume; ie2++)
        for(int ie4 = 0; ie4< nume; ie4++)
        {
            prop11->SetDeltaDM(ie4,ie2);
            dm2t11=prop11->PropagateDM(times2);
   
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
	    transition.AssignLevels(0,ie1+numg,ig1);
            //ids[0] = edipsdata2D[ig1][ie1+numg];

            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

        for(int ig2 = 0; ig2< numg; ig2++)
        {
	    transition.AssignLevels(1,ie4+numg,ig2);
            //ids[1] = edipsdata2D[ie4+numg][ig2];
        for(int ie5 = 0; ie5< nume; ie5++)
        {
	    transition.AssignLevels(2,ie5+numg,if1+nume_numg);
            //ids[2] = edipsdata2D[ie5+numg][if1+nume_numg];
        for(int ie6 = 0; ie6< nume; ie6++)
        for(int if2 = 0; if2< numf; if2++)
        {

	    transition.AssignLevels(3,ie6+numg,if2+nume_numg);
            //ids[3] = edipsdata2D[ie6+numg][if2+nume_numg];

            // result
	    complexv damplitude = transition.GetAveragedAmplitude();
//            double damplitude = obj_avor.rot_av_dip_4(ies,ids);
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3n.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *conj(dm1t10.data3D[it1][ie2][ig2])
                    *dm2t11.data3D[it2][ie5][ie3]
                    *dm3t21.data3D[it3][if2][ie6];
            }
        }}}}}}
    }
    
    // all K1 diagrams have been calculated
    
    /////////////////
    // K2 diagrams
    
    if(coherentGSBK2 || transportGSBK2)
    { 
        if(prop00 == 0)
	{
            prop00 = new Type(*this);
	    	prop00->flagNonsecular = nonsecular;
		prop00->flagMarkovian=markovian;
		prop00->calcR->flagLindblad = lindblad;
		prop00->SetBlock(0);
	}

        if(prop10 == 0)
	{
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
		prop10->flagMarkovian=markovian;
		prop10->calcR->flagLindblad = lindblad;
		prop10->SetBlock(10);
	}
    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

     // coherent GSB K2
    //     ig6 |  | ig6
    //        -------
    //     ie3 |  | ig6
    //     ie2 |  | ig5
    //        -------
    //     ig3 |  | ig5
    //     ig4 |  | ig2
    //        -------
    //     ie4 |  | ig2
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
	    transition.AssignLevels(0,ie1+numg,ig1);
            //ids[0] = edipsdata2D[ig1][ie1+numg];

            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

            //if(!markovian){
            //    prop10->PropagateDM(times2);
            //    dm2auxt10=prop10->GetHistory();
            //}
            
        for(int ig2 = 0; ig2< numg; ig2++)
        for(int ig4 = 0; ig4< numg; ig4++)
        {
            prop00->SetDeltaDM(ig4,ig2);
            dm2t00=prop00->PropagateDM(times2);
   
        for(int ig5 = 0; ig5< numg; ig5++)
        for(int ie2 = 0; ie2< nume; ie2++)
        {
            prop10->SetDeltaDM(ie2,ig5);
            //if(!markovian){
            //    prop10->SetHistory(dm2auxt10);
            //}
            dm3t10=prop10->PropagateDM(times3);

        for(int ie4 = 0; ie4< nume; ie4++)
        {
	    transition.AssignLevels(1,ie4+numg,ig4);
            //ids[1] = edipsdata2D[ie4+numg][ig4];
        for(int ig3 = 0; ig3< numg; ig3++)
        {
	    transition.AssignLevels(2,ie2+numg,ig3);
            //ids[2] = edipsdata2D[ig3][ie2+numg];
        for(int ie3 = 0; ie3< nume; ie3++)
        for(int ig6 = 0; ig6< numg; ig6++)
        {

	    transition.AssignLevels(3,ie3+numg,ig6);
            //ids[3] = edipsdata2D[ie3+numg][ig6];

            // result
	    complexv damplitude = transition.GetAveragedAmplitude();
//            double damplitude = obj_avor.rot_av_dip_4(ies,ids);
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3d.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *dm1t10.data3D[it1][ie4][ig2]
                    *dm2t00.data3D[it2][ig3][ig5]
                    *dm3t10.data3D[it3][ie3][ig6];
            }
        }}}}}}
    }
    
    if(coherentESEK2 || transportESEK2)
    { 

        if(prop10 == 0)
	{
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
		prop10->flagMarkovian=markovian;
		prop10->calcR->flagLindblad = lindblad;
		prop10->SetBlock(10);
	}
        if(prop11 == 0)
	{
            prop11 = new Type(*this);
	    	prop11->flagNonsecular = nonsecular;
		prop11->flagMarkovian=markovian;
		prop11->calcR->flagLindblad = lindblad;
		prop11->SetBlock(11);
	}

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

    // coherent ESE k2
    //     ig4 |  | ig4
    //        -------
    //     ie6 |  | ig4
    //     ie5 |  | ig3
    //        -------
    //     ie5 |  | ie3
    //     ie2 |  | ie4
    //        -------
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
    
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
	    transition.AssignLevels(0,ie1+numg,ig1);
            //ids[0] = edipsdata2D[ig1][ie1+numg];

            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

            //if(!markovian){
            //    prop10->PropagateDM(times2);
            //    dm2auxt10=prop10->GetHistory();
            //}

        for(int ie2 = 0; ie2< nume; ie2++)
        for(int ie4 = 0; ie4< nume; ie4++)
        {
            prop11->SetDeltaDM(ie2,ie4);
            dm2t11=prop11->PropagateDM(times2);
   
        for(int ig3 = 0; ig3< numg; ig3++)
        for(int ie5 = 0; ie5< nume; ie5++)
        {
            prop10->SetDeltaDM(ie5,ig3);
            //if(!markovian){
            //    prop10->SetHistory(dm2auxt10);
            //}
            dm3t10=prop10->PropagateDM(times3);

        for(int ig2 = 0; ig2< numg; ig2++)
        {
	    transition.AssignLevels(1,ie4+numg,ig2);
            //ids[1] = edipsdata2D[ie4+numg][ig2];
        for(int ie3 = 0; ie3< nume; ie3++)
        {
	    transition.AssignLevels(2,ie3+numg,ig3);
            //ids[2] = edipsdata2D[ie3+numg][ig3];
        for(int ie6 = 0; ie6< nume; ie6++)
        for(int ig4 = 0; ig4< numg; ig4++)
        {

	    transition.AssignLevels(3,ie6+numg,ig4);
            //ids[3] = edipsdata2D[ie6+numg][ig4];

            // result
	    complexv damplitude = transition.GetAveragedAmplitude();
//            double damplitude = obj_avor.rot_av_dip_4(ies,ids);
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3d.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *dm1t10.data3D[it1][ie2][ig2]
                    *dm2t11.data3D[it2][ie5][ie3]
                    *dm3t10.data3D[it3][ie6][ig4];
            }
        }}}}}}
    }
      
    if(coherentESAK2 || transportESAK2)
    { 


        if(prop10 == 0)
	{
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
		prop10->flagMarkovian=markovian;
		prop10->calcR->flagLindblad = lindblad;
		prop10->SetBlock(10);
	}
        if(prop11 == 0)
	{
            prop11 = new Type(*this);
	    	prop11->flagNonsecular = nonsecular;
		prop11->flagMarkovian=markovian;
		prop11->calcR->flagLindblad = lindblad;
		prop11->SetBlock(11);
	}
        if(prop21 == 0)
	{
            prop21 = new Type(*this);
	    	prop21->flagNonsecular = nonsecular;
		prop21->flagMarkovian=markovian;
		prop21->calcR->flagLindblad = lindblad;
		prop21->SetBlock(21);
	}

    	ios[0] = xomega1;
    	ios[1] = -xomega2;
    	ios[2] = xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

    // coherent ESA k2
    //     ie6 |  | ie6
    //        -------
    //     if2 |  | ie6
    //     if1 |  | ie3
    //        -------
    //     ie5 |  | ie3
    //     ie4 |  | ie2
    //        -------
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
        for(int if1 = 0; if1< numf; if1++)
        for(int ie3 = 0; ie3< nume; ie3++)
        {
            prop21->SetDeltaDM(if1,ie3);
            dm3t21=prop21->PropagateDM(times3);

        for(int ie2 = 0; ie2< nume; ie2++)
        for(int ie4 = 0; ie4< nume; ie4++)
        {
            prop11->SetDeltaDM(ie4,ie2);
            dm2t11=prop11->PropagateDM(times2);
   
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
	    transition.AssignLevels(0,ie1+numg,ig1);
            //ids[0] = edipsdata2D[ig1][ie1+numg];

            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

        for(int ig2 = 0; ig2< numg; ig2++)
        {
	    transition.AssignLevels(1,ie2+numg,ig2);
            //ids[1] = edipsdata2D[ie2+numg][ig2];
        for(int ie5 = 0; ie5< nume; ie5++)
        {
	    transition.AssignLevels(2,ie5+numg,if1+nume_numg);
            //ids[2] = edipsdata2D[ie5+numg][if1+nume_numg];
        for(int ie6 = 0; ie6< nume; ie6++)
        for(int if2 = 0; if2< numf; if2++)
        {

	    transition.AssignLevels(3,ie6+numg,if2+nume_numg);
            //ids[3] = edipsdata2D[ie6+numg][if2+nume_numg];

            // result
	    complexv damplitude = transition.GetAveragedAmplitude();
//            double damplitude = obj_avor.rot_av_dip_4(ies,ids);
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3n.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *dm1t10.data3D[it1][ie2][ig2]
                    *dm2t11.data3D[it2][ie5][ie3]
                    *dm3t21.data3D[it3][if2][ie6];
            }
        }}}}}}
    }
    
    
   
    
    
    //////// now KIII
    
    if(coherentES1K3)
    { 

        if(prop10 == 0)
	{
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
		prop10->flagMarkovian=markovian;
		prop10->calcR->flagLindblad = lindblad;
		prop10->SetBlock(10);
	}
        if(prop20 == 0)
	{
            prop20 = new Type(*this);
	    	prop20->flagNonsecular = nonsecular;
		prop20->flagMarkovian=markovian;
		prop20->calcR->flagLindblad = lindblad;
		prop20->SetBlock(20);
	}
    	ios[0] = xomega1;
    	ios[1] = xomega2;
    	ios[2] = -xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

    // coherent ESA1 k3
    //     ig4 |  | ig4
    //        -------
    //     ie4 |  | ig4
    //     ie3 |  | ig3
    //        -------
    //     if2 |  | ig3
    //     if1 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
	    transition.AssignLevels(0,ie1+numg,ig1);
            //ids[0] = edipsdata2D[ig1][ie1+numg];

            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

            //if(!markovian){
            //    prop10->PropagateDM(times2);
            //    dm2auxt10=prop10->GetHistory();
            //}
            
        for(int if1 = 0; if1< numf; if1++)
        for(int ig2 = 0; ig2< numg; ig2++)
        {
            prop20->SetDeltaDM(if1,ig2);
            dm2t20=prop20->PropagateDM(times2);
   
        for(int ie3 = 0; ie3< nume; ie3++)
        for(int ig3 = 0; ig3< numg; ig3++)
        {
            prop10->SetDeltaDM(ie3,ig3);
            //if(!markovian){
            //    prop10->SetHistory(dm2auxt10);
            //}
            dm3t10=prop10->PropagateDM(times3);

        for(int ie2 = 0; ie2< nume; ie2++)
        {
	    transition.AssignLevels(1,ie2+numg,if1+nume_numg);
            //ids[1] = edipsdata2D[ie2+numg][if1+nume_numg];
        for(int if2 = 0; if2< numf; if2++)
        {
	    transition.AssignLevels(2,ie3+numg,if2+nume_numg);
            //ids[2] = edipsdata2D[ie3+numg][if2+nume_numg];
        for(int ie4 = 0; ie4< nume; ie4++)
        for(int ig4 = 0; ig4< numg; ig4++)
        {

	    transition.AssignLevels(3,ie4+numg,ig4);
            //ids[3] = edipsdata2D[ie4+numg][ig4];

            // result
	    complexv damplitude = transition.GetAveragedAmplitude();
//            double damplitude = obj_avor.rot_av_dip_4(ies,ids);
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3d.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *dm1t10.data3D[it1][ie2][ig2]
                    *dm2t20.data3D[it2][if2][ig3]
                    *dm3t10.data3D[it3][ie4][ig4];
            }
        }}}}}}
    }
    
    if(coherentES2K3)
    { 

        if(prop10 == 0)
	{
            prop10 = new Type(*this);
	    	prop10->flagNonsecular = nonsecular;
		prop10->flagMarkovian=markovian;
		prop10->calcR->flagLindblad = lindblad;
		prop10->SetBlock(10);
	}
        if(prop20 == 0)
	{
            prop20 = new Type(*this);
	    	prop20->flagNonsecular = nonsecular;
		prop20->flagMarkovian=markovian;
		prop20->calcR->flagLindblad = lindblad;
		prop20->SetBlock(20);
	}
        if(prop21 == 0)
	{
            prop21 = new Type(*this);
	    	prop21->flagNonsecular = nonsecular;
		prop21->flagMarkovian=markovian;
		prop21->calcR->flagLindblad = lindblad;
		prop21->SetBlock(21);
	}
    	ios[0] = xomega1;
    	ios[1] = xomega2;
    	ios[2] = -xomega3;
    	ios[3] = -xomega4;
	transition.PopulateFields(ies,iks,ios);

    // coherent ESA2 k3
    //     ie4 |  | ie4
    //        -------
    //     if3 |  | ie4
    //     if2 |  | ie3
    //        -------
    //     if2 |  | ig3
    //     if1 |  | ig2
    //        -------
    //     ie2 |  | ig2
    //     ie1 |  | ig1
    //        -------
    //     ig1 |  | ig1
        for(int if2 = 0; if2< numf; if2++)
        for(int ie3 = 0; ie3< nume; ie3++)
        {
            prop21->SetDeltaDM(if2,ie3);
            dm3t21=prop21->PropagateDM(times3);

        for(int if1 = 0; if1< numf; if1++)
        for(int ig2 = 0; ig2< numg; ig2++)
        {
            prop20->SetDeltaDM(if1,ig2);
            dm2t20=prop20->PropagateDM(times2);
   
        for(int ig1 = 0; ig1< numg; ig1++)
        for(int ie1 = 0; ie1< nume; ie1++)
        {
	    transition.AssignLevels(0,ie1+numg,ig1);
            //ids[0] = edipsdata2D[ig1][ie1+numg];

            prop10->SetDeltaDM(ie1,ig1);
            dm1t10=prop10->PropagateDM(times1);

        for(int ie2 = 0; ie2< nume; ie2++)
        {
	    transition.AssignLevels(1,ie2+numg,if1+nume_numg);
            //ids[1] = edipsdata2D[ie2+numg][if1+nume_numg];
        for(int ig3 = 0; ig3< numg; ig3++)
        {
	    transition.AssignLevels(2,ie3+numg,ig3);
            //ids[2] = edipsdata2D[ie3+numg][ig3];
        for(int ie4 = 0; ie4< nume; ie4++)
        for(int if3 = 0; if3< numf; if3++)
        {

	    transition.AssignLevels(3,ie4+numg,if3+nume_numg);
            //ids[3] = edipsdata2D[ie4+numg][if3+nume_numg];

            // result
            //double damplitude = obj_avor.rot_av_dip_4(ies,ids);
	    complexv damplitude = transition.GetAveragedAmplitude();
            // time loops
            for(int it1=0; it1<time1N; it1++)
            for(int it2=0; it2<time2N; it2++)
            for(int it3=0; it3<time3N; it3++)
            {
            dat3n.data3D[it3][it2][it1] += cnni*damplitude
                    *grPopsdata1D[ig1]
                    *dm1t10.data3D[it1][ie2][ig2]
                    *dm2t20.data3D[it2][if2][ig3]
                    *dm3t21.data3D[it3][if3][ie4];
            }
        }}}}}}
    } 


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
}



