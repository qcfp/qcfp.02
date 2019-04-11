#include"feinman2sideddiagram3.hpp"
#include"../averages_orient/averages_orient.hpp"
//#include"../cuda_calculate/cuda_calculate.hpp"

feinman2sideddiagram3::feinman2sideddiagram3()
{
    cout<<"Error: you should not use feinman2sideddiagram3()\n";
}
feinman2sideddiagram3::feinman2sideddiagram3(complexv* ifreq, dvector3d* ids, dvector3d* ies,int& averaging)
{
 
    // energy for the three intervals (imaginary part is for the dephasing)
	freq = ifreq;

	// amplitudes of transitions (dipoles)
	//ds = ids;

	// vectors of fields (dipoles)
	//es = ies;

        averages_orient obj_avor(averaging);
        damplitude = obj_avor.rot_av_dip_4(ies,ids);
        
    ongrF = 0;
    ongij = 0;
    
    technique = -1;
    states = 0;
    //for(int is = 0; is< 14; is ++)
    //{
    //    states[is]=-1;
    //}

    feinman2sideddiagram3smallparameter = -1;
    flagCheck = 0;

}
feinman2sideddiagram3::feinman2sideddiagram3(complexv* ifreq, complexv factor)
{
    
    // energy for the three intervals (imaginary part is for the dephasing)
    freq = ifreq;
    
    damplitude = factor;
    
    ongrF = 0;
    ongij = 0;
    
    technique = -1;
    states = 0;
    //for(int is = 0; is< 14; is ++)
    //{
    //    states[is]=-1;
    //}
    
    
    feinman2sideddiagram3smallparameter = -1;
    flagCheck = 0;
    
}

feinman2sideddiagram3::~feinman2sideddiagram3()
{
        if(ongij)
            delete c4p_obj;// = new correlation4points(gij,gijAmplitude,taud,tauc,taub,taua);
        ongij = 0;
        ongrF = 0;

}


void feinman2sideddiagram3::propagate(storage<complexv>& res,double* itT, int* itN)
{
    int& time3N = itN[2];
    int& time2N = itN[1];
    int& time1N = itN[0];

    double& time3i = itT[5];
    double& time2i = itT[4];
    double& time1i = itT[3];
    double time3s = (itT[2]-itT[5])/time3N;
    double time2s = (itT[1]-itT[4])/time2N;
    double time1s = (itT[0]-itT[3])/time1N;
    double tim3;
    double tim2;
    double tim1;

    complexv*** dat = res.data3D;
//---------------------------------------------------------------------------------
//--- GPU CALCULATION -------------------------------------------------------------
//---------------------------------------------------------------------------------
/*    gRe = gf->DirectAccessRE();
    gIm = gf->DirectAccessIM();
    gRdat = gRe->DirectAccess();
    gIdat = gIm->DirectAccess();
    gRD = gRdat->data1D;
    gID = gIdat->data1D;
    int gfN = gRdat->GetSize();
    double gfstep = gRe->GetStep();
    double gfmin = gRe->GetXI();
    double gfmax = gRe->GetXF();
    double grmin, grmax, grstep;
    int grn;
    if (ongrF)
    {
	grdat = grF->DirectAccess();
	grd = grdat->data1D;
	grn = grdat->GetSize();
        grstep = grF->GetStep();
        grmin = grF->GetXI();
        grmax = grF->GetXF();
    }
    data_real = new double **[time3N];
    data_imag = new double **[time3N];
    freq_real = new double [3];
    freq_imag = new double [3];
    for (int i = 0; i < 3; i++)
    {
	freq_real[i] = freq[i].real();
	freq_imag[i] = freq[i].imag();
    }
    for (int i = 0; i < time3N; i++)
    {
        data_real[i] = new double *[time2N];
        data_imag[i] = new double *[time2N];
        for (int j = 0; j < time2N; j++)
        {
            data_real[i][j] = new double [time1N];
            data_imag[i][j] = new double [time1N];
        }
    }
//    for (int itime1 = 0; itime1 < time1N; itime1++)
//	for (int itime2 = 0; itime2 < time2N; itime2++)
//	    for (int itime3 = 0; itime3 < time3N; itime3++)
//	    {
//		data_real[itime3][itime2][itime1] = dat[itime3][itime2][itime1].real();
//		data_imag[itime3][itime2][itime1] = dat[itime3][itime2][itime1].imag();
//            }
//    cout << "Technique: " << technique << "\n";
    int cusuccess = cuda_calculate(time3N, time2N, time1N, time3s, time2s, time1s, data_real, data_imag, ongrF, freq_real, freq_imag, ongij, technique, gRD, gID, gfN, grd, grn, gfamp, nstates, gfmin, gfmax, gfstep, grmin, grmax, grstep, damplitude, inumO, taud, tauc, taub, taua);
    for (int i = 0; i < time3N; i++)
	for (int j = 0; j < time2N; j++)
	    for (int k = 0; k < time1N; k++)
	    {
		complexv temporary(data_real[i][j][k], data_imag[i][j][k]);
		dat[i][j][k] += temporary;
//		cout << dat[i][j][k] << "\n";
            }
    for (int i = 0; i < time3N; i++)
    {
	for (int j = 0; j < time2N; j++)
	{
	    delete [] data_real[i][j];
	    delete [] data_imag[i][j];
        }
	delete [] data_real[i];
	delete [] data_imag[i];
    }
    delete [] data_real;
    delete [] data_imag;
    delete [] freq_real;
    delete [] freq_imag;*/

//---------------------------------------------------------------------------------
//------ GPU CALCULATION END ------------------------------------------------------
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//------ CPU CALCULATION ----------------------------------------------------------
//---------------------------------------------------------------------------------
        for(int itime2 = 0; itime2<time2N; itime2++)// loop over time 2
        {
            tim2 =  time2i+time2s * itime2;
    for(int itime3 = 0; itime3<time3N; itime3++)// loop over time 3
    {
        tim3 =  time3i+time3s * itime3;
            for(int itime1 = 0; itime1<time1N; itime1++)// loop over time 1
            {
                tim1 =  time1i+time1s * itime1;
                
                if(flagCheck==1)
                {
                    //cout<<"flag ckecked\n";
                    if(time2N == 1)
                        return;
                }
          //      cout<<itime3<<" "<<itime2<<" "<<itime1<<" "<<dat[itime3][itime2][itime1]<<"\n";

                dat[itime3][itime2][itime1] += calculate(tim3,tim2,tim1);
        
        //cout<<itime3<<" "<<itime2<<" "<<itime1<<" "<<dat[itime3][itime2][itime1]<<"\n";
            }
        }
    }

//---------------------------------------------------------------------------------
//----------- CPU CALCULATION END -------------------------------------------------
//---------------------------------------------------------------------------------
}


        
complexv feinman2sideddiagram3::calculate(double& tim3, double& tim2, double& tim1)
{
	complexv amp(0.0,0.0);
        //cout<<amp<<"  test \n";
	if(ongrF==1)
	{
            double intensity = grF->Get(tim2);
            if(intensity<feinman2sideddiagram3smallparameter)
            {
                flagCheck = 1;
                return 0.0;
            }
            // here is incoherent population diagram
	        amp = cnni*(freq[2]*tim3+freq[0]*tim1);
       
                if(ongij)
                {
                    if(technique == 9) // kii
                    {
                        amp+= c4p_obj->GetIncoherent(tim3,tim2,tim1);
                    }
                    else if(technique == 10) //ki
                    {
                        amp+= conj( c4p_obj->GetIncoherent(tim3,tim2,tim1) );
                    }
                }
                
                

                amp = exp(amp);
                amp = amp*intensity;

	}
	if(ongrF==2)
	{
            double intensity = propF->GetAt(tim2);
            //cout<<"intensity: "<<"\n";
            if(intensity<feinman2sideddiagram3smallparameter)
            {
                flagCheck = 1;
                return 0.0;
            }

            // here is incoherent population diagram
	        amp = cnni*(freq[2]*tim3+freq[0]*tim1);

                if(ongij)
                {
                    if(technique == 9) // kii
                    {
                        amp+= c4p_obj->GetIncoherent(tim3,tim2,tim1);
                    }
                    else if(technique == 10) //ki
                    {
                        amp+= conj( c4p_obj->GetIncoherent(tim3,tim2,tim1) );
                        //amp+=  c4p_obj->GetIncoherent(tim3,tim2,tim1);
                    }
                }



                amp = exp(amp);
                amp = amp*intensity;

	}
	else
	{
            // here is the coherent diagram
            amp = cnni*(freq[2]*tim3+freq[1]*tim2+freq[0]*tim1);
            
//            cout<<"coherent A "<<freq[2]<<" "<<freq[1]<<" "<<freq[0]<<"\n";    

            
            if(ongij)
            {
                double tau4,tau3,tau2,tau1;
                if(technique ==1) // llll
                {
                   tau4 = tim1+tim2+tim3;
                   tau3 = tim1+tim2;
                   tau2 = tim1;
                   tau1 = 0;
                }
                else if(technique ==2) // lllr
                {
//                    cout<<"tim3="<<tim3<<"\n";
//                    cout<<"tim2="<<tim2<<"\n";
//                    cout<<"tim1="<<tim1<<"\n";
//                    
                   tau4 = 0;
                   tau3 = tim1+tim2+tim3;
                   tau2 = tim1+tim2;
                   tau1 = tim1;
                }
                else if(technique ==3) // llrl
                {
                   tau4 = tim1;
                   tau3 = tim1+tim2+tim3;
                   tau2 = tim1+tim2;
                   tau1 = 0;
                }
                else if(technique ==4) // llrr
                {
                   tau4 = 0;
                   tau3 = tim1;
                   tau2 = tim1+tim2+tim3;
                   tau1 = tim1+tim2;
                }
                else if(technique ==5) // lrll
                {
                   tau4 = tim1+tim2;
                   tau3 = tim1+tim2+tim3;
                   tau2 = tim1;
                   tau1 = 0;
                }
                else if(technique ==6) // lrlr
                {
                   tau4 = 0;
                   tau3 = tim1+tim2;
                   tau2 = tim1+tim2+tim3;
                   tau1 = tim1;
                }
                else if(technique ==7) // lrrl
                {
                   tau4 = tim1;
                   tau3 = tim1+tim2;
                   tau2 = tim1+tim2+tim3;
                   tau1 = 0;
                }
                else if(technique ==8) // lrrrr
                {
                   tau4 = 0;
                   tau3 = tim1;
                   tau2 = tim1+tim2;
                   tau1 = tim1+tim2+tim3;
                }

           
                    amp+= c4p_obj->Get(tau4,tau3,tau2,tau1);

             }
             
            
                amp = exp(amp);
        }
            
            
        //cout<<"coherent B "<<amp<<"\n";    
        
        
        
        amp = amp* damplitude;

        
        return amp;

}



void feinman2sideddiagram3::assignLineshape(int& itechnique,asymptoticLF_complexv* igij,double*** igijAmplitude,int& inumOsc,int* istates, int& ninput)
{
//    numOsc = inumOsc;

    
    //for(int is = 0; is< 14; is ++)
    //{
        states=istates;
    //}
    
    technique = itechnique;

    
    // coherent diagrams
               if(technique ==1) // llll
                {
                   taud = states[3];
                   tauc = states[2];
                   taub = states[1];
                   taua = states[0];
                }
                else if(technique ==2) // lllr
                {
                   taud = states[5];
                   tauc = states[3];
                   taub = states[2];
                   taua = states[0];
                }
                else if(technique ==3) // llrl
                {
                   taud = states[5];
                   tauc = states[3];
                   taub = states[1];
                   taua = states[0];
                }
                else if(technique ==4) // llrr
                {
                   taud = states[7];
                   tauc = states[5];
                   taub = states[3];
                   taua = states[0];
                }
                else if(technique ==5) // lrll
                {
                   taud = states[5];
                   tauc = states[2];
                   taub = states[1];
                   taua = states[0];
                }
                else if(technique ==6) // lrlr
                {
                   taud = states[6];
                   tauc = states[5];
                   taub = states[2];
                   taua = states[0];
                }
                else if(technique ==7) // lrrl
                {
                   taud = states[6];
                   tauc = states[5];
                   taub = states[1];
                   taua = states[0];
                }
                else if(technique ==8) // lrrrr
                {
                   taud = states[7];
                   tauc = states[6];
                   taub = states[5];
                   taua = states[0];
                }
    
        // incoherent diagrams
                else if(technique ==9) // kii 
                {
                   taud = states[9];
                   tauc = states[5];
                   taub = states[3];
                   taua = states[1];
                }
                if(technique ==10) // ki
                {
                   taud = states[9];
                   tauc = states[3];
                   taub = states[5];
                   taua = states[7];
                }

    
    c4p_obj = new correlation4points(igij,igijAmplitude,inumOsc,taud,tauc,taub,taua);
    gf = igij;
    inumO = inumOsc;
    gfamp = igijAmplitude;
    ongij = 1;
    nstates = ninput;
    
    
    
}

void feinman2sideddiagram3::assignPropagator(asymptoticLF<double>* iGij)
{
    // this propagator is the incoherent Pauli master propagator 
    // from state c into state j
    grF = iGij;
    ongrF = 1;
}
void feinman2sideddiagram3::assignPropagator(propagatorM* iGij)
{
    // this propagator is the incoherent Pauli master propagator
    // from state c into state j
    propF = iGij;
    ongrF = 2;
}

//void feinman2sideddiagram3::multiplier(double& ival)
//{
//    damplitude *= ival;
//}
