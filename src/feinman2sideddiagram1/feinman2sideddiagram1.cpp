#include"feinman2sideddiagram1.hpp"
#include"../averages_orient/averages_orient.hpp"
	
feinman2sideddiagram1::feinman2sideddiagram1(complexv& ifreq, interaction& itr)
{
 
    // energy (imaginary part is for the dephasing)
	freq = ifreq;


    damplitude = itr.GetAveragedAmplitude();

    gij = 0;
}
feinman2sideddiagram1::feinman2sideddiagram1(complexv& ifreq, dvector3d* ids, dvector3d* ies,int& averaging)
{
 
    // energy (imaginary part is for the dephasing)
	freq = ifreq;

	// amplitudes of transitions (dipoles)
	ds = ids;

	// vectors of fields (dipoles)
	es = ies;
        
    averages_orient obj_avor(averaging);
    damplitude = obj_avor.rot_av_dip_2(es,ds);

    gij = 0;
}
feinman2sideddiagram1::feinman2sideddiagram1(complexv& ifreq, double& idamplitude)
{
    
    // energy (imaginary part is for the dephasing)
    freq = ifreq;
    
    damplitude = idamplitude;
    
    gij = 0;

}

feinman2sideddiagram1::~feinman2sideddiagram1()
{
}


void feinman2sideddiagram1::propagate(interpolationF<complexv>& res)
{
    int timeN = res.GetN();
    double timeT = res.GetXF();
    double times = timeT/timeN;
    double tim;
    complexv* dat = res.DirectAccessD();
   
    for(int itime = 0; itime<timeN; itime++)// loop over time
	{
		tim =  times * itime;
                
        dat[itime]=calculate(tim);
                //cout<<dat[itime]<<"\n";
    }
}

void feinman2sideddiagram1::propagateC(interpolationF<complexv>& res)
{
    int timeN = res.GetN();
    double timeT = res.GetXF();
    double times = timeT/timeN;
    double tim;
    complexv* dat = res.DirectAccessD();
    
    for(int itime = 0; itime<timeN; itime++)// loop over time
    {
        tim =  times * itime;
        
        dat[itime]=(calculate(tim)).conjugate();
        //cout<<dat[itime]<<"\n";
    }
}



complexv feinman2sideddiagram1::calculate(double& tim)
{
        complexv amp = -coni*freq*tim;

        if(gij != 0)
        {
        	for(int io = 0; io<numOsc; io ++)
        		amp -= (gijAmplitude[io]+gijAmplitude2[io]-gijAmplitude3[io]-gijAmplitude4[io])*gij[io].Get(tim);
        }
        amp = exp(amp);
        amp *= damplitude;
        
        return  coni*amp;
}

void feinman2sideddiagram1::assignLineshape(asymptoticLF_complexv* igij,double* igijAmplitude,int& inumOsc)
{
//    gij = igij;
//    gijAmplitude = igijAmplitude;
//    numOsc = inumOsc;
    
    cout<<"ERROR: this function is valid only for a single ground state!!!\n";
}

void feinman2sideddiagram1::assignLineshape(asymptoticLF_complexv* igij,double* igijA1,double* igijA2,double* igijA3,double* igijA4,int& inumOsc)
{
    gij = igij;
    gijAmplitude = igijA1;
    gijAmplitude2 = igijA2;
    gijAmplitude3 = igijA3;
    gijAmplitude4 = igijA4;
    numOsc = inumOsc;
}
