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

#include"calculator_3rd_secular_cumulant.hpp"



#ifdef MPIPROCESSING
#include<mpi.h>
#endif

calculator_3rd_secular_cumulant::calculator_3rd_secular_cumulant()
:communicator_3rd()
{
    //Setup();
    eigensys = true;
    outputstring = "# QCFP calculator_3rd_secular_cumulant class signal\n";
    
    storedat1.SetDimension(3);
    storedat2.SetDimension(3);
    storedat3.SetDimension(3);
    
    transport = 1;
    speedupsmallness = -1;
    noComplexLifetimes = 0;
    nonmarkovian = 1;
    
}    
calculator_3rd_secular_cumulant::calculator_3rd_secular_cumulant(string& ifname)
:communicator_3rd()
{
    //Setup();
    eigensys = true;
    outputstring = "# QCFP calculator_3rd_secular_cumulant class signal\n";
    
    storedat1.SetDimension(3);
    storedat2.SetDimension(3);
    storedat3.SetDimension(3);
    
    transport = 1;
    speedupsmallness = -1;
    noComplexLifetimes = 0;
    nonmarkovian = 1;
    
    
    
    // reading all options from the input file
    toolsIO tio;
   	tio.ReadWholeFile(ifname,inputfile,1);

    ifstream istr(ifname.c_str());
    if(!istr.is_open())
    {
        cout<<"Error: input file cannot be read\nquitting.\n";
        return;
    }
    tio.StreamSkipTrailers(istr);
    string tstr; getline(istr,tstr);
    ReadSystem(istr);
    ReadBath(istr);
    ReadExperiment3rd(istr);
	// specific methods for this calculator
	storage<double> arr(2);
	storage<int> iri(2);
 
	if(tio.LookUpAndReadSquareMatrix<int>(
		"MethodSpeedUpSmallness:",
		"Reading the speed-up parameter\n", 
		1,1, iri,inputfile.str()))
	{
        	speedupsmallness = iri.data2D[0][0];
	}
	
    std::size_t found;

    found = experiment_complexity.find("NoComplexLifetimes");
    if(found!=std::string::npos)
	noComplexLifetimes = 1;
    
    found = experiment_complexity.find("Markovian");
    if(found!=std::string::npos)
	nonmarkovian = 0;
    
    found = experiment_complexity.find("NonMarkovian");
    if(found!=std::string::npos)
	nonmarkovian = 1;
    
    // done reading. Next will do computations
   
}


void calculator_3rd_secular_cumulant::LaunchBasic()
{
    // here all information must be specified so that the object knows how to continue
    if(simType == 1)
    {
        cout<<"Warning: TWW is only for 2QC\n";
        signal=LaunchTWW();
    }
    else if(simType == 2)
    {
        cout<<"Warning: WTW is only for RP and NRP\n";
        signal=LaunchWTW();
    }
    else if(simType == 3)
    {
        cout<<"Warning: WWT is only for  2QC\n";
        signal=LaunchWWT();
    }
    else if(simType == 4)
    {
        cout<<"Warning: TTT is for full time domain response function\n";
        LaunchTTT();
    }
    else
    {
        cout<<"Error in calculator_3rd_secular_cumulant: simulation type is improper\n";
    }
}
interpolationF2d<complexv> calculator_3rd_secular_cumulant::LaunchTWW(double iifre2, double iffre2, double iifre3, double iffre3, int inump, double itf1)
{
    ifre3 = iifre3;
    ffre3 = iffre3;
    ifre2 = iifre2;
    ffre2 = iffre2;
    nump = inump;
    tf1 = itf1;
    return LaunchTWW();
}

interpolationF2d<complexv> calculator_3rd_secular_cumulant::LaunchWTW(double iifre1, double iffre1, double iifre3, double iffre3, int inump, double itf2)
{
    ifre1 = iifre1;
    ffre1 = iffre1;
    ifre3 = iifre3;
    ffre3 = iffre3;
    nump = inump;
    tf2 = itf2;
    return LaunchWTW();
}
interpolationF2d<complexv> calculator_3rd_secular_cumulant::LaunchWWT(double iifre1, double iffre1, double iifre2, double iffre2, int inump, double itf3)
{
    ifre1 = iifre1;
    ffre1 = iffre1;
    ifre2 = iifre2;
    ffre2 = iffre2;
    nump = inump;
    tf3 = itf3;
    return LaunchWWT();
}


interpolationF2d<complexv> calculator_3rd_secular_cumulant::LaunchWTW()
// depends on ifre1, ffre1, ifre3, ffre3, nump, tf2
{
    cout<<"Launching WTW\n";
    constants cst;
    
 
    if( (coherentGSBK1 == coherentGSBK2)&&(coherentESEK1 == coherentESEK2)&&(coherentESAK1 == coherentESAK2) && (transportGSBK1 == transportGSBK2)&&(transportESEK1 == transportESEK2)&&(transportESAK1 == transportESAK2))
        {
            // for empty calculation
            if( (coherentGSBK1 == 0)&&(coherentESEK1 == 0)&&(coherentESAK1 == 0) && (transportGSBK1 == 0)&&(transportESEK1 == 0)&&(transportESAK1 == 0))
            {
                return interpolationF2d<complexv>();
            }
        

            // this is pp simulation
            // this is complete without optimization as if two separate runs

            // saving pattern
            int lcb1 =coherentGSBK1;
            int lce1 =coherentESEK1;
            int lca1 =coherentESAK1;
            int lcb2 =coherentGSBK2;
            int lce2 =coherentESEK2;
            int lca2 =coherentESAK2;
            int ltb1 =transportGSBK1;
            int lte1 =transportESEK1;
            int lta1 =transportESAK1;
            int ltb2 =transportGSBK2;
            int lte2 =transportESEK2;
            int lta2 =transportESAK2;

            // these must be positive
            double lfi1 = ifre1;
            double lff1 = ffre1;

            // first run
            coherentGSBK2 = 0;
            transportGSBK2 = 0;
            coherentESEK2 = 0;
            transportESEK2 = 0;
            coherentESAK2 = 0;
            transportESAK2 = 0;
            
            ifre1 = -lfi1;
            ffre1 = -lff1;
            
            interpolationF2d<complexv> sigk1 =LaunchWTW();
            
            // first run
            coherentGSBK1 = 0;
            transportGSBK1 = 0;
            coherentESEK1 = 0;
            transportESEK1 = 0;
            coherentESAK1 = 0;
            transportESAK1 = 0;

            coherentGSBK2 = lcb2;
            transportGSBK2 = ltb2;
            coherentESEK2 = lce2;
            transportESEK2 = lte2;
            coherentESAK2 = lca2;
            transportESAK2 = lta2;

            ifre1 = lfi1;
            ffre1 = lff1;
            
            interpolationF2d<complexv> sigk2 =LaunchWTW();

        
        
        
        
            // restoring
            coherentGSBK1 = lcb1;
            transportGSBK1 = ltb1;
            coherentESEK1 = lce1;
            transportESEK1 = lte1;
            coherentESAK1 = lca1;
            transportESAK1 = lta1;
            
            
            
            // combining the result
            complexv** id1 = sigk1.DirectAccessD();
            complexv** id2 = sigk2.DirectAccessD();
            
            for(int i3=0; i3<nump; i3++)
            for(int i1=0; i1<nump; i1++)
                {
                    id2[i1][i3] += id1[i1][i3];
                }
            
            return sigk2;

        
        }
    else
    {
        
        double cfre1 = 0.5*(ifre1+ffre1);
        double freS1 = (ffre1-ifre1)/nump;
        double freF1 = ffre1-ifre1;
        
        double cfre3 = 0.5*(ifre3+ffre3);
        double freS3 = (ffre3-ifre3)/nump;
        double freF3 = ffre3-ifre3;
        
        
        // making time variables
        int timeN = nump/2;
        double timeS3 = fabs(cst.pi2/freF3);
        double timeF3 = fabs(timeS3*timeN);
        double timeS1 = fabs(cst.pi2/freF1);
        double timeF1 = fabs(timeS1*timeN);
        
        double timeDP[6];
        int timeIP[3];
        // final values
        timeDP[0]=timeF1;
        timeDP[1]=tf2;
        timeDP[2]=timeF3;
        // initial values
        timeDP[3]=0;
        timeDP[4]=tf2;
        timeDP[5]=0;
        // numbers
        timeIP[0]=timeN;
        timeIP[1]=1;
        timeIP[2]=timeN;
        
        //cout<<"timeDP[4]=tf2; "<<timeDP[4]<<"\n";
        
        storage<complexv> dat3d(3);
        dat3d.Allocate(timeN,1,timeN);
        // this one is for positive contributions
        
        storage<complexv> dat3n(3);
        dat3n.Allocate(timeN,1,timeN);
        // this one is for negative contributions
        
        
        // making time domain simulations
        LaunchGeneric(dat3d,dat3n,timeDP,timeIP);
        
        // now one has to make FFT
        
        // moving oscillations to zero frequency:
        
        
        for(int it1=0; it1<timeN; it1++)
            for(int it3=0; it3<timeN; it3++)
            {
                double t1 = timeS1*it1;
                double t3 = timeS3*it3;
                
                complexv phase = exp(coni*(cfre1*t1+cfre3*t3)-naturallinewidth*(t1+t3));
                dat3d.data3D[it3][0][it1] = phase*dat3d.data3D[it3][0][it1];
                dat3n.data3D[it3][0][it1] = phase*dat3n.data3D[it3][0][it1];
            }
        
        
        //making 2D response functions
        complexv** idat;
        interpolationF2d<complexv> asymres(0,timeF1/timeN,timeN,  0,timeF3/timeN,timeN);
        idat = asymres.DirectAccessD();
        
        for(int i1=0; i1<timeN; i1++)
            for(int i3=0; i3<timeN; i3++)
            {
                idat[i1][i3] = dat3d.data3D[i3][0][i1]-dat3n.data3D[i3][0][i1];
            }
        dat3d.Delete();
        dat3n.Delete();
        
        // scaling zero points due to Fourier transformation
        for(int id=0; id<timeN; id++)
        {
            idat[0][id] /= 2.0;
            idat[id][0] /= 2.0;
        }
        
        //    cout<<asymres.Get(0)<<"\n";
        //    cout<<asymres.Get(0.1)<<"\n";
        
        
        // doing Fourier
        interpolationF2d<complexv> asymresFFTi(0,timeS1,nump,  0,timeS3,nump);
        interpolationF2d<complexv> asymresFFTf(0,freS1,nump,  0,freS3,nump);
        asymresFFTi *= complexv(0,0);
        asymresFFTi.DirectAssign(asymres);
        
        //    cout<<asymresFFTi.Get(0)<<"\n";
        //    cout<<asymresFFTi.Get(0.1)<<"\n";
        
        
        toolsFFT fft;
        asymresFFTf = fft.executeN(asymresFFTi);
        
        
        //    cout<<asymresFFTf.Get(0)<<"\n";
        //    cout<<asymresFFTf.Get(0.1)<<"\n";
        
        // shifting frequencies
        fft.SwapSides(asymresFFTf);
        
        if(freS1<0)
        {
            // switching sides to make the same as in kII
            complexv** id1 = asymresFFTf.DirectAccessD();
            interpolationF2d<complexv> tmpf = asymresFFTf;
            complexv** id2 = tmpf.DirectAccessD();
            
            for(int i3=0; i3<nump; i3++)
                for(int i1=0; i1<nump; i1++)
                {
                    id1[i1][i3] = id2[nump-i1-1][i3];
                }
            
        }
        asymresFFTf.UpdateAxis(ifre1,freS1, ifre3,freS3);

        // output
        return asymresFFTf;
        
    }
}


interpolationF2d<complexv> calculator_3rd_secular_cumulant::LaunchTWW()
// depends on ifre2, ffre2, ifre3, ffre3, nump, tf1
{
    cout<<"Launching TWW\n";
    constants cst;
    
    double cfre3 = 0.5*(ifre3+ffre3);
    double cfre2 = 0.5*(ifre2+ffre2);
    
    double freS3 = (ffre3-ifre3)/nump;
    double freF3 = ffre3-ifre3;
    double freS2 = (ffre2-ifre2)/nump;
    double freF2 = ffre2-ifre2;
    
    
    // making time variables
    int timeN = nump/2;
    double timeS3 = fabs(cst.pi2/freF3);
    double timeF3 = fabs(timeS3*timeN);
    double timeS2 = fabs(cst.pi2/freF2);
    double timeF2 = fabs(timeS2*timeN);
    
    double timeDP[6];
    int timeIP[3];
    // numbers
    timeIP[0]=1;
    timeIP[1]=timeN;
    timeIP[2]=timeN;
    // final values
    timeDP[0]=tf1;
    timeDP[1]=timeF2;
    timeDP[2]=timeF3;
    // initial values
    timeDP[3]=tf1;
    timeDP[4]=0;
    timeDP[5]=0;
    
    
    storage<complexv> dat3d(3);
    dat3d.Allocate(timeN,timeN,1);
    // this one is for positive contributions
    
    storage<complexv> dat3n(3);
    dat3n.Allocate(timeN,timeN,1);
    // this one is for negative contributions
    
    
    // making time domain simulations
    LaunchGeneric(dat3d,dat3n,timeDP,timeIP);    
    
    // now one has to make FFT
    
    // moving oscillations:
    
    for(int it2=0; it2<timeN; it2++)
        for(int it3=0; it3<timeN; it3++)
        {
            double t2 = timeS2*it2;
            double t3 = timeS3*it3;
            
            complexv phase = exp(coni*(cfre2*t2+cfre3*t3)-naturallinewidth*(t2+t3));
            dat3d.data3D[it3][it2][0] *= phase;
            dat3n.data3D[it3][it2][0] *= phase;
        }
    
    
    //making 2D response functions
    complexv** idat;
    interpolationF2d<complexv> asymres(0,timeF2/timeN,timeN,  0,timeF3/timeN,timeN);
    idat = asymres.DirectAccessD();
    
    for(int i2=0; i2<timeN; i2++)
        for(int i3=0; i3<timeN; i3++)
        {
            idat[i2][i3] = dat3d.data3D[i3][i2][0]-dat3n.data3D[i3][i2][0];
        }
    dat3d.Delete();
    dat3n.Delete();
    
    // scaling zero points due to Fourier transformation
    for(int id=0; id<timeN; id++)
    {
        idat[0][id] /= 2.0;
        idat[id][0] /= 2.0;
    }
    
    //    cout<<asymres.Get(0)<<"\n";
    //    cout<<asymres.Get(0.1)<<"\n";
    
    
    // doing Fourier
    interpolationF2d<complexv> asymresFFTi(0,timeS2,nump,  0,timeS3,nump);
    interpolationF2d<complexv> asymresFFTf(0,freS2,nump,  0,freS3,nump);
    asymresFFTi *= complexv(0,0);
    asymresFFTi.DirectAssign(asymres);
    
    //    cout<<asymresFFTi.Get(0)<<"\n";
    //    cout<<asymresFFTi.Get(0.1)<<"\n";
    
    
    toolsFFT fft;
    asymresFFTf = fft.executeN(asymresFFTi);
    
    
    //    cout<<asymresFFTf.Get(0)<<"\n";
    //    cout<<asymresFFTf.Get(0.1)<<"\n";
    
    // shifting frequencies
    fft.SwapSides(asymresFFTf);
    asymresFFTf.UpdateAxis(ifre2,freS2, ifre3,freS3);
    
    
    //    cout<<asymresFFTf.Get(0)<<"\n";
    //    cout<<asymresFFTf.Get(0.1)<<"\n";
    
    // output 
    return asymresFFTf;
}


interpolationF2d<complexv> calculator_3rd_secular_cumulant::LaunchWWT()
// depends on ifre1, ffre1, ifre2, ffre2, nump, tf3
{
    cout<<"Launching WWT\n";

    constants cst;
    
    double cfre1 = 0.5*(ifre1+ffre1);
    double cfre2 = 0.5*(ifre2+ffre2);
    
    double freS1 = (ffre1-ifre1)/nump;
    double freF1 = ffre1-ifre1;
    double freS2 = (ffre2-ifre2)/nump;
    double freF2 = ffre2-ifre2;
    
    
    // making time variables
    int timeN = nump/2;
    double timeS1 = fabs(cst.pi2/freF1);
    double timeF1 = fabs(timeS1*timeN);
    double timeS2 = fabs(cst.pi2/freF2);
    double timeF2 = fabs(timeS2*timeN);
    
    double timeDP[6];
    int timeIP[3];
    
    // numbers
    timeIP[0]=timeN;
    timeIP[1]=timeN;
    timeIP[2]=1;
    
    // final values
    timeDP[0]=timeF1;
    timeDP[1]=timeF2;
    timeDP[2]=tf3;

    // initial values
    timeDP[3]=0;
    timeDP[4]=0;
    timeDP[5]=tf3;

    
    
    storage<complexv> dat3d(3);
    dat3d.Allocate(1,timeN,timeN);
    // this one is for positive contributions
    
    storage<complexv> dat3n(3);
    dat3n.Allocate(1,timeN,timeN);
    // this one is for negative contributions
    
    
    // making time domain simulations
    LaunchGeneric(dat3d,dat3n,timeDP,timeIP);    
    
    // now one has to make FFT
    
    // moving oscillations:
    
    for(int it2=0; it2<timeN; it2++)
        for(int it1=0; it1<timeN; it1++)
        {
            double t2 = timeS2*it2;
            double t1 = timeS1*it1;
            
            complexv phase = exp(coni*(cfre2*t2+cfre1*t1)-naturallinewidth*(t1+t2));
            dat3d.data3D[0][it2][it1] *= phase;
            dat3n.data3D[0][it2][it1] *= phase;
        }
    
    
    //making 2D response functions
    complexv** idat;
    interpolationF2d<complexv> asymres(0,timeF1/timeN,timeN,  0,timeF2/timeN,timeN);
    idat = asymres.DirectAccessD();
    
    for(int i2=0; i2<timeN; i2++)
        for(int i1=0; i1<timeN; i1++)
        {
            idat[i1][i2] = dat3d.data3D[0][i2][i1]-dat3n.data3D[0][i2][i1];
        }
    dat3d.Delete();
    dat3n.Delete();
    
    // scaling zero points due to Fourier transformation
    for(int id=0; id<timeN; id++)
    {
        idat[0][id] /= 2.0;
        idat[id][0] /= 2.0;
    }
    
    //    cout<<asymres.Get(0)<<"\n";
    //    cout<<asymres.Get(0.1)<<"\n";
    
    
    // doing Fourier
    interpolationF2d<complexv> asymresFFTi(0,timeS1,nump,  0,timeS2,nump);
    interpolationF2d<complexv> asymresFFTf(0,freS1,nump,  0,freS2,nump);
    asymresFFTi *= complexv(0,0);
    asymresFFTi.DirectAssign(asymres);
    
    //    cout<<asymresFFTi.Get(0)<<"\n";
    //    cout<<asymresFFTi.Get(0.1)<<"\n";
    
    
    toolsFFT fft;
    asymresFFTf = fft.executeN(asymresFFTi);
    
    
    //    cout<<asymresFFTf.Get(0)<<"\n";
    //    cout<<asymresFFTf.Get(0.1)<<"\n";
    
    // shifting frequencies
    fft.SwapSides(asymresFFTf);
    asymresFFTf.UpdateAxis(ifre1,freS1, ifre2,freS2);
    
    
    //    cout<<asymresFFTf.Get(0)<<"\n";
    //    cout<<asymresFFTf.Get(0.1)<<"\n";
    
    // output 
    return asymresFFTf;
}


void calculator_3rd_secular_cumulant::LaunchTTT()
// depends on ffre1, ffre2, ti1, ti2, ti3, tf1, tf2, tf3, nump
{
    cout<<"Launching TTT\n";

    constants cst;
    
    // these are used as central frequencies
    double cfre1 = 0;
    double cfre2 = 0;
    double cfre3 = 0;
    
    
    // making time variables
    int timeN = nump;
    
    double timeS3 = (tf3-ti3)/nump;
    double timeS2 = (tf2-ti2)/nump;
    double timeS1 = (tf1-ti1)/nump;
    
    double timeDP[6];
    int timeIP[3];
    
    // numbers
    timeIP[0]=timeN;
    timeIP[1]=timeN;
    timeIP[2]=timeN;
    
    // final values
    timeDP[0]=tf1;
    timeDP[1]=tf2;
    timeDP[2]=tf3;

    // initial values
    timeDP[3]=ti1;
    timeDP[4]=ti2;
    timeDP[5]=ti3;
    
    storedat1.Allocate(nump,nump,nump);
    storedat2.Allocate(nump,nump,nump);
    storedat3.Allocate(nump,nump,nump);

    
    storage<complexv> dat3d(3);
    dat3d.Allocate(timeN,timeN,timeN);
    // this one is for positive contributions
    
    storage<complexv> dat3n(3);
    dat3n.Allocate(timeN,timeN,timeN);
    // this one is for negative contributions

            // saving pattern
            int lcb1 =coherentGSBK1;
            int lce1 =coherentESEK1;
            int lca1 =coherentESAK1;
            int lcb2 =coherentGSBK2;
            int lce2 =coherentESEK2;
            int lca2 =coherentESAK2;
            int ltb1 =transportGSBK1;
            int lte1 =transportESEK1;
            int lta1 =transportESAK1;
            int ltb2 =transportGSBK2;
            int lte2 =transportESEK2;
            int lta2 =transportESAK2;
            int lcf1 = coherentES1K3;
            int lcf2 = coherentES2K3;

            
            
    // running rephasing
            coherentGSBK1 = lcb1;
            coherentESEK1 = lce1;
            coherentESAK1 = lca1;
            coherentGSBK2 = 0;
            coherentESEK2 = 0;
            coherentESAK2 = 0;
            transportGSBK1 = ltb1;
            transportESEK1 = lte1;
            transportESAK1 = lta1;
            transportGSBK2 = 0;
            transportESEK2 = 0;
            transportESAK2 = 0;
            coherentES1K3 = 0;
            coherentES2K3 = 0;
    
     cfre1 = -ifre1;
     cfre2 = ifre2-ifre1;
     cfre3 = ifre3+ifre2-ifre1;
    
    
    // making time domain simulations
    LaunchGeneric(dat3d,dat3n,timeDP,timeIP);    
    
    // moving oscillations:
    for(int it3=0; it3<timeN; it3++)
      for(int it2=0; it2<timeN; it2++)
        for(int it1=0; it1<timeN; it1++)
        {
            double t3 = timeS3*it3;
            double t2 = timeS2*it2;
            double t1 = timeS1*it1;
            
            complexv phase = exp(coni*(cfre3*t3+cfre2*t2+cfre1*t1)-naturallinewidth*(t1+t2+t3));
            storedat1.data3D[it3][it2][it1] = (dat3d.data3D[it3][it2][it1]-dat3n.data3D[it3][it2][it1])*phase;
            
            dat3d.data3D[it3][it2][it1] = 0.0;
            dat3n.data3D[it3][it2][it1] = 0.0;
        }

    
    // running non-rephasing
            coherentGSBK1 = 0;
            coherentESEK1 = 0;
            coherentESAK1 = 0;
            coherentGSBK2 = lcb2;
            coherentESEK2 = lce2;
            coherentESAK2 = lca2;
            transportGSBK1 = 0;
            transportESEK1 = 0;
            transportESAK1 = 0;
            transportGSBK2 = ltb2;
            transportESEK2 = lte2;
            transportESAK2 = lta2;
            coherentES1K3 = 0;
            coherentES2K3 = 0;
    
     cfre1 = ifre1;
     cfre2 = -ifre2+ifre1;
     cfre3 = ifre3-ifre2+ifre1;
    
    
    // making time domain simulations
    LaunchGeneric(dat3d,dat3n,timeDP,timeIP);    
    
    // moving oscillations:
    for(int it3=0; it3<timeN; it3++)
      for(int it2=0; it2<timeN; it2++)
        for(int it1=0; it1<timeN; it1++)
        {
            double t3 = timeS3*it3;
            double t2 = timeS2*it2;
            double t1 = timeS1*it1;
            
            complexv phase = exp(coni*(cfre3*t3+cfre2*t2+cfre1*t1)-naturallinewidth*(t1+t2+t3));
            storedat2.data3D[it3][it2][it1] = (dat3d.data3D[it3][it2][it1]-dat3n.data3D[it3][it2][it1])*phase;
            
            dat3d.data3D[it3][it2][it1] = 0.0;
            dat3n.data3D[it3][it2][it1] = 0.0;
        }

    
      // running double quantum coherence
            coherentGSBK1 = 0;
            coherentESEK1 = 0;
            coherentESAK1 = 0;
            coherentGSBK2 = 0;
            coherentESEK2 = 0;
            coherentESAK2 = 0;
            transportGSBK1 = 0;
            transportESEK1 = 0;
            transportESAK1 = 0;
            transportGSBK2 = 0;
            transportESEK2 = 0;
            transportESAK2 = 0;
            coherentES1K3 = lcf1;
            coherentES2K3 = lcf2;
    
     cfre1 = ifre1;
     cfre2 = ifre2+ifre1;
     cfre3 = -ifre3+ifre2+ifre1;
    
    
    // making time domain simulations
    LaunchGeneric(dat3d,dat3n,timeDP,timeIP);    
    
    // moving oscillations:
    for(int it3=0; it3<timeN; it3++)
      for(int it2=0; it2<timeN; it2++)
        for(int it1=0; it1<timeN; it1++)
        {
            double t3 = timeS3*it3;
            double t2 = timeS2*it2;
            double t1 = timeS1*it1;
            
            complexv phase = exp(coni*(cfre3*t3+cfre2*t2+cfre1*t1)-naturallinewidth*(t1+t2+t3));
            storedat3.data3D[it3][it2][it1] = (dat3d.data3D[it3][it2][it1]-dat3n.data3D[it3][it2][it1])*phase;
            
            dat3d.data3D[it3][it2][it1] = 0.0;
            dat3n.data3D[it3][it2][it1] = 0.0;
        }
  
    
                // restoring  all:
    
            coherentGSBK1 = lcb1;
            coherentESEK1 = lce1;
            coherentESAK1 = lca1;
            coherentGSBK2 = lcb2;
            coherentESEK2 = lce2;
            coherentESAK2 = lca2;
            transportGSBK1 = ltb1;
            transportESEK1 = lte1;
            transportESAK1 = lta1;
            transportGSBK2 = ltb2;
            transportESEK2 = lte2;
            transportESAK2 = lta2;
            coherentES1K3 = lcf1;
            coherentES2K3 = lcf2;


    

    dat3d.Delete();
    dat3n.Delete();

    // output 
    return;
}





void calculator_3rd_secular_cumulant::PublishSpecial4(ofstream& ofs)
{
    // this is only for rephasing and non-rephasing 2D
        ofs<<outputstring;
        ofs.precision(12);

        double tim1;
        double tim2;
        double tim3;
        double dt1 = (tf1-ti1)/nump;
        double dt2 = (tf2-ti2)/nump;
        double dt3 = (tf3-ti3)/nump;
        
        
        ofs<<ifre1<<"\t"<<ifre2<<"\t"<<ifre3<<"\t"<<nump<<"\n";
        
        for(int ind1 = 0; ind1<nump; ind1++)
            for(int ind2 = 0; ind2<nump; ind2++)
            for(int ind3 = 0; ind3<nump; ind3++)
            {
            tim1 = ti1+dt1*ind1;
            tim2 = ti2+dt2*ind2;
            tim3 = ti3+dt3*ind3;

            ofs<<tim3<<"\t"<<tim2<<"\t"<<tim1<<"\t";

            complexv result = storedat1.data3D[ind3][ind2][ind1];
            ofs<<result.real()<<"\t";
            ofs<<result.imag()<<"\t";

            result = storedat2.data3D[ind3][ind2][ind1];
            ofs<<result.real()<<"\t";
            ofs<<result.imag()<<"\t";
            
            result = storedat3.data3D[ind3][ind2][ind1];
            ofs<<result.real()<<"\t";
            ofs<<result.imag()<<"\n";
            
        }

}


