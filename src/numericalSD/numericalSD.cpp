// everything is performed in terms of the accessory bath functions a1,a2 and a3
// at the moment that algorithm is most advanced
#include"../storage/storage.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../interpolationF/interpolationF.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../complexv/complexv.hpp"
#include"../toolsFFT/toolsFFT.hpp"
#include"../constants/constants.hpp"

#include"numericalSD.hpp"


#include<string>
#include<iomanip>
using std::setprecision;

numericalSD::numericalSD()
{
    constants cc;

    const0_695 = 1.0;
    const5305 = 1.0;
    const3000 = cc.pi2*const5305;

    temperature = 0;

    spectral_density = 0;
	accessory_a1P = 0;
	accessory_a2P = 0;
	accessory_a3P = 0;
    lineshape = 0;

    globalgamma = 0.0;
    epsilonomega = 0.0;
    flagaccessories = 0;
}

numericalSD::numericalSD( double dtemperature)
{
    constants cc;

    const0_695 = 1.0;
    const5305 = 1.0;
    const3000 = cc.pi2*const5305;

    temperature = const0_695*dtemperature;

    spectral_density = 0;
    accessory_a1P = 0;
    accessory_a2P = 0;
    accessory_a3P = 0;
    lineshape = 0;


    globalgamma = 0.0;
    epsilonomega = 0.0;
    flagaccessories = 0;

    // reading the sepctral density from file
    //ReadSpectralDensities(file_name);

}

numericalSD::numericalSD( double dtemperature,asymptoticLF<double>& ifun)
{
    //cout<<"inside SD\n";
    constants cc;

    const0_695 = 1.0;
    const5305 = 1.0;
    const3000 = cc.pi2*const5305;

    temperature = const0_695*dtemperature;

    spectral_density = 0;
	accessory_a1P = 0;
	accessory_a2P = 0;
	accessory_a3P = 0;
    lineshape = 0;

    globalgamma = 0.0;
    epsilonomega = 0.0;
    flagaccessories = 0;


   Allocate_spectral_density();
    *spectral_density = ifun;
    spectral_density->causality = 3; // odd function

	int numv = 2*spectral_density->GetN();
	double xmax = 2*spectral_density->GetXF(); // cm-1
	epsilonomega = xmax/numv/50; // cm-1



	// reading the sepctral density from file
	//ReadSpectralDensities(file_name);

}

void numericalSD::updateUnits(string istring)
{
    constants cc;

    // reconverting temperature
	temperature = temperature/const0_695;

    if(istring == "default")
    {

        const0_695 = 1.0;
        const5305  = 1.0;
        const3000 = cc.pi2*const5305;
    }
    else if(istring == "cm-fs-K")
    {

        const0_695 = cc.bkTc;
        const5305  = cc.c_cmfs_i;
        const3000 = cc.pi2*const5305;
    }


    // converting temperature
	temperature = temperature*const0_695;
}




numericalSD::~numericalSD()
{
	Delete_spectral_density();
    //cout<<"1 stuch here?\n";
	DeleteMBARFunction();
    //cout<<"2 stuch here?\n";
	DeleteAccessoryFunction();
    //cout<<"3 stuch here?\n";
	DeleteLineshape();
    //cout<<"4 stuch here?\n";
}



void numericalSD::Allocate_spectral_density()
{
	// here we allocate empty spectral densities
	// we use it only for positive time
	// because for negative time is essentially symmetric

	if(spectral_density == NULL)
		spectral_density = new asymptoticLF<double>;

	else
		cout<<"Error: error with spectral density\n";
}


void numericalSD::AllocateAccessoryFunction()
{
	// call to this function is OK when reading SD is done
	if(spectral_density != NULL && accessory_a1P == NULL)
	{
		accessory_a1P = new asymptoticLF_complexv;
		accessory_a2P = new asymptoticLF_complexv;
		accessory_a3P = new asymptoticLF_complexv;

	}
	else
	{
		cout<<"Error: error with spectral density\n";
	}
}
void numericalSD::AllocateLineshape()
{
	// call to this function is OK when reading SD is done
	if(spectral_density != NULL && lineshape == NULL)
	{
		lineshape = new asymptoticLF_complexv;
        lineshape->SetCausality(2,3); // causal function
	}
	else
	{
		cout<<"Error: error with spectral density\n";
	}
}
void numericalSD::AllocateMBARFunction()
// this is half-FFT of the correlation function in the positive side
{
//no need for allocation
}

void numericalSD::DeleteAccessoryFunction()
{
// these are complex functions

	if(accessory_a1P != NULL)
	{
		delete accessory_a1P;// = new AsymptoticLinearFunction_complexv;
		delete accessory_a2P;// = new AsymptoticLinearFunction_complexv;
		delete accessory_a3P;// = new AsymptoticLinearFunction_complexv;
	}
	accessory_a1P = NULL;
	accessory_a2P = NULL;
	accessory_a3P = NULL;
}
void numericalSD::DeleteLineshape()
{
// these are complex functions

	if(lineshape != NULL)
	{
		delete lineshape;
        }
	lineshape = NULL;
}

void numericalSD::DeleteMBARFunction()
{
// no need to delete
}

void numericalSD::Delete_spectral_density()
{
	if(spectral_density!=NULL)
		delete spectral_density;
	spectral_density = NULL;
}


void numericalSD::SetSD(double* idat,double dx,int nump)
{
   	if(spectral_density == NULL)
		Allocate_spectral_density();

        storage<double> tdat(1);
        tdat.Allocate(nump);
        idat = tdat.FlipBar(idat);

	*spectral_density = asymptoticLF<double>(tdat,0,dx);//MakeFromDataset(idat,dx,nump);
        idat = tdat.FlipBar(idat);
        tdat.Delete();
        spectral_density->a=0;
        spectral_density->b=0;
        spectral_density->c=0;
        spectral_density->d=0;
        spectral_density->causality = 3;

    	int numv = 2*spectral_density->GetN();
    	double xmax = 2*spectral_density->GetXF(); // cm-1
    	epsilonomega = xmax/numv/50; // cm-1



}
void numericalSD::SetSD(asymptoticLF<double>& ifun)
{
    //cout<<"inside SD\n";
   if(spectral_density == NULL)
       Allocate_spectral_density();
    *spectral_density = ifun;
    spectral_density->causality = 3;

	int numv = 2*spectral_density->GetN();
	double xmax = 2*spectral_density->GetXF(); // cm-1
	epsilonomega = xmax/numv/50; // cm-1

}
void numericalSD::ReadSD(string inp_file_densities)
{
	// reads spectral density from file in numerical format

	if(spectral_density == NULL)
		Allocate_spectral_density();


	if (inp_file_densities != "" )
	{
			ifstream inp_file_str(inp_file_densities.c_str());

			if(!inp_file_str.is_open())
				cout<<"Error: bath spectral densities file has not been openned\n";
			else
			{
				cout<<"Reading bath spectral densities\n";

                asymptoticLF<double> spdens;
                interpolationF<double> spdensRd;
                spdensRd.ReadF(&inp_file_str);
                spdens = spdensRd;
                spdens.causality = 3;

                SetSD(spdens);
			}

	}

	int numv = 2*spectral_density->GetN();
	double xmax = 2*spectral_density->GetXF(); // cm-1
	epsilonomega = xmax/numv/50; // cm-1


}

void numericalSD::ReadCfun(string inp_file)
{
	// reads spectral density from file in numerical format

	if(spectral_density == NULL)
		Allocate_spectral_density();


	if (inp_file != "" )
	{
			ifstream inp_file_str(inp_file.c_str());

			if(!inp_file_str.is_open())
				cout<<"Error: the classical correlation function file has not been openned\n";
			else
			{
				cout<<"Reading bath classical correlation function\n";

                asymptoticLF<double> spdens;
                interpolationF<double> spdensRd;
                spdensRd.ReadF(&inp_file_str);

                // now I have an interpolation function - which is the correlation function
                int nump = spdensRd.GetN();
                double step = spdensRd.GetStep(); // this guy is in outer units

                // 1-st step : making the C' function; converting units
                interpolationF<complexv> ifun(0.0,step/const5305,2*nump);
                complexv* hifun = ifun.DirectAccessD();
                double* hspdensRd = spdensRd.DirectAccessD();
                hifun[0]=hspdensRd[0];
                for(int ind=1;ind<nump;ind++)
                {
                    hifun[ind]=hspdensRd[ind];
                    hifun[2*nump-ind]=hifun[ind];
                }

                toolsFFT tfft;
                interpolationF<complexv> c1 = tfft.executeN(ifun);
                hifun = c1.DirectAccessD();
                double wstep=c1.GetStep();

                // making spectral density
                for(int ind=0;ind<nump;ind++)
                {
                    hspdensRd[ind]=(hifun[ind]).real()*step/const5305*tanh(ind*wstep/temperature/2.0);
                }
                spdensRd.UpdateAxis(0,wstep);

                spdens = spdensRd;
                spdens.causality = 3;

                SetSD(spdens);
			}

	}

	int numv = 2*spectral_density->GetN();
	double xmax = 2*spectral_density->GetXF(); // cm-1
	epsilonomega = xmax/numv/50; // cm-1


}

double numericalSD::GetSD(double dfreq)
// cfreq in wavenumbers (not angular)
// returns function in wavenumber units
{

    if(spectral_density == 0)
    {
        cout<<"ERROR in numericalSD: spectral density is empty\n";
        return 0;
    }

    // this function accepts just real frequencies
	//double sign = 1.0;
	//double val;
	//if(dfreq<0.0)
	//{
	//	dfreq = -dfreq;
	//	sign  = -sign;
	//}
    //    if(dfreq > (spectral_density->GetXF()-spectral_density->GetXF()/spectral_density->GetN()) )
    //        return 0;

	// dfreq is always positive
    //val = spectral_density->Get(dfreq);
	//return sign*val;

    return spectral_density->Get(dfreq);
}

complexv numericalSD::GetMf(double omega)
{


    if(accessory_a1P == NULL)
        ConstructAccessoryFunction();
	// numerical form of M function must be created



    // this function takes some default values and calculates spectral density integrals
    double valuer;


	int numv = 2*spectral_density->GetN();
	double xmax = 2*spectral_density->GetXF(); // cm-1
	double deltav = xmax/numv; // cm-1

	double&kT = temperature;


    // real part
	if(fabs(omega)<0.001*deltav)
		valuer = kT*GetSD(deltav)/deltav;
	else if(fabs(100*omega)<kT)
		valuer = 0.5*GetSD(omega)*(2*kT/omega + 1);
	else
        valuer = 0.5*GetSD(omega)*(1.0/tanh(0.5*omega/kT)+1.0);


   //imaginary part


	return complexv(valuer, M_BAR_im.Get(omega));


}

complexv numericalSD::GetMfn(double omega)
{
    return GetMf(-omega);
}
complexv numericalSD::GetMfp(double omega)
{
    return GetMf(omega);
}



complexv numericalSD::GetGfD(double t, int n)
{
	// this function returns derivatives of the lineshape function
	// input time is in outer units
	// output is
	//    n=0: dimensionless
	//    n=1: cm-1
	//    n=2: cm-2
	// n-th (1 and 2) derivative

	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();
	if(lineshape == NULL)
		CreateLineshape();

    t=t/const5305;

	double tempt=fabs(t);

	complexv result;

	if(n==0) return GetGf(t);

	else if( n == 1 )
	{
		result =  coni*accessory_a2P->Get(tempt)
			-coni*accessory_a2P->Get(0);
		if (t>=0)
			return result;
		else
			return -conj(result);
		}
	else if( n == 2 )
	{
            // this is the correlation function
		result = accessory_a3P->Get(tempt);
		if (t>=0)
			return result;
		else
			return conj(result);
	}
	else
	{
		cout<<"error: only up to second derivative is available\n";
		return 0.0;
	}
}

// returns lineshape function
complexv numericalSD::GetGf(double t)
{
    t=t/const5305;

    if(accessory_a1P == NULL)
		ConstructAccessoryFunction();
	if(lineshape == NULL)
		CreateLineshape();

//	double tempt=fabs(t);


    // causality is included
    return lineshape->Get(t);
}

// returns lineshape function from accessory functions (totally relies on FFT; may be wrong at long times)
complexv numericalSD::GetGfA(double t)
{
    t=t/const5305;

    if(accessory_a1P == NULL)
		ConstructAccessoryFunction();

	double tempt=fabs(t);


	// using accessory functions if available (this is one of the most accurate and fastest numerical results
	complexv gfun = accessory_a1P->Get(0)
			-accessory_a1P->Get(tempt)
			-coni*accessory_a2P->Get(0)*tempt/const5305;


    if(t<0.0) return conj(gfun);
	else return gfun;
}




complexv numericalSD::GetGfD1Inf()
{

//    cout<<"Fix here to use M-funcion\n";
    //cout<<"mfun: "<< GetMf(0)<<"\n";
    return GetMf(0);


    // temporary functions
    //asymptoticLF<double>* rep = lineshape->DirectAccessRE();
    //double vr = rep->a;
    //rep = lineshape->DirectAccessIM();
    //double vi = rep->a;
    //cout<<"afun: "<< complexv(vr,vi)<<"\n";
    //return complexv(vr,vi);

}




////////////////////////////////
//  internal functions



void  numericalSD::CreateLineshape()
{
    if(lineshape == NULL)
        AllocateLineshape();
    else
        return;


    if(accessory_a1P == NULL)
        ConstructAccessoryFunction();

    complexv* gfun = 0;
    int numt=0;
    double deltat=0;




    if(flagaccessories < 2) // can be 0 or 1
    {
        cout<<"Creating lineshape function: accessory functions\n";
        // use accessory
        //	// using accessory functions if available (this is one of the most accurate and fastest numerical results
        //	complexv gfun = accessory_a1P->GetValueApproach("0")
        //			-accessory_a1P->GetValueApproach(tempt)
        //			-coni*accessory_a2P->GetValueApproach("0")*tempt/const5305;

        //	// using direct lineshape
        //	complexv gfun = lineshape->Get(tempt);
        //
        //        if(t<0.0) return conj(gfun);
        //	else return gfun;


        // update
        numt = accessory_a3P->GetN();
        double tmax = accessory_a3P->GetXF();
        //cout<<"tma: "<<numt<<" "<<tmax<<"\n";

        deltat = tmax/numt;

        complexv lam = -coni*accessory_a2P->Get(0)*deltat;
        gfun = new complexv[numt+4];

        for(int id = 0; id<numt+4; id++)
            gfun[id] = accessory_a1P->Get(0)-accessory_a1P->Get(id*deltat)+lam*((double)id);



    }
    else
    {

        // accessory is now used only for the correlation function

        cout<<"Creating lineshape function: direct integration\n";

        int numv = spectral_density->GetN();

        // Using smart integration of the time correlation function

        int factor = 1;

        numt = factor*numv;
        double vmax = factor*spectral_density->GetXF(); // cm-1

        double deltav = vmax/numt; // cm-1
        deltat = const3000/(2*vmax); // cm  = 2*pi*1

        double dt2 = deltat*deltat;

        gfun = new complexv[numt+4];
        complexv* Cort = new complexv[numt+4];
        complexv  cval;

        // assigning gfun:
        for(int ind=0;ind<numt+4;ind++)
            Cort[ind]=GetCf(ind*deltat);


        // beginning the iterative integration
        gfun[0] = 0.0;

        for( int indi = 1; indi < numt+4; indi++)
        {
            //  calculating Zn
            cval =  0.5*Cort[0];

            for(int indj = 1; indj < indi-1; indj++) //(* first integral positive *)
                cval +=  Cort[indj];


            cval +=  Cort[indi-1]*0.875;

            cval +=  Cort[indi]*0.125;

            gfun[indi] = cval*dt2+gfun[indi-1];
        }
        //(* that was only positive time part - negative will be conjugate *)

        delete[] Cort;


        //        cout<<numt<<" test\n";
    }

    // assigning the lineshape function
    storage<complexv> tgfun(1);

    // adding ZPL linewidth
    if(gfun[numt+1].real()<gfun[numt].real())
    	gfun[numt+1] = complexv(gfun[numt].real(),gfun[numt+1].imag());

    for(int ind=0;ind<numt+4;ind++)
    	gfun[ind]=gfun[ind] + (ind*deltat)*globalgamma;

    complexv derivative = (gfun[numt+1]-gfun[numt])/deltat;

    tgfun.Allocate(numt);
    tgfun.SetBar(gfun,numt);
    //    *lineshape = asymptoticLF_complexv(tgfun,0,deltat);

    //cout<<"gfu: "<<numt<<" "<<deltat<<"\n";


    *lineshape = asymptoticLF_complexv(tgfun,0,deltat,
                                       complexv(0.0,0.0),
                                       gfun[numt],derivative);


    lineshape->SetCausality(2,3);
    tgfun.Delete();
    delete[] gfun;

    //GetGfD1Inf();

}


void  numericalSD::ConstructAccessoryFunction()
{
	// Here I calculate very accurately. The Lineshape function should be
	// always represented through the accessory functions.

//	if(spectral_density == NULL)
//		ReadSpectralDensities();
	if(accessory_a1P == NULL)
	    	AllocateAccessoryFunction();
	else
		return;

	cout<<"Creating accessory bath functions\n";

	// This function can be used only for numerical spectral density representation
	// When lineshape function derivatives are being used
	// here I calculate accessory temporal functions based on the
	// quantum correlation function in frequency domain
	// dummySD[ind]
	// that correlation function is organized as for the Fourier transform




	// Using FFT to create time domain functions
        // a1 and a2 will have small time step
        // function a3 will have timestep consistent with fft

    // 6 times longer FFT - from experience
    // no smaller than 2
    // int fact = 1; // 6;

    // basic length is from the spectral density
    int numv = spectral_density->GetN();
	double vmax = spectral_density->GetXF(); // cm-1
	double deltav = vmax/numv; // cm-1

    // including negative, we must have at least 2*
    int numt = 2*numv;
    double deltat = const3000/(2*vmax); // fs  = 2*pi*5305

	double* dummySD = new double[numt];// double length
	double* dummyPP = new double[numt];// double length

    // the following function creates proper array ordering for FFT (positive and negative)
    CreateFDCorrelationSD(dummySD,dummyPP,deltav,numt);
    double ksi = dummyPP[1]/deltav;

	complexv* funa1 = new complexv[numt];
	complexv* funa2 = new complexv[numt];
	complexv* funa3 = new complexv[numt];
    complexv* times = new complexv[numt];

    // arrays funax include the negative part of the time as well

        //cout<<funa1[0]<<"\n";
        //cout<<funa1[1]<<"\n";
        //cout<<funa1[2]<<"\n";


    for(int ind = 0; ind<numt; ind++)
        times[ind] = deltat*ind;



    if(flagaccessories == 1)
    {
        // doing complete run
        GenerateAccessoryArrays(dummySD, dummyPP, funa1, funa2, funa3, numt, deltav);
    }

    else
    {
        // doing normal run for A3 function
        GenerateAccessoryArrays(dummySD, dummyPP, 0, 0, funa3, numt, deltav);

        // A1 and A2 functions are obtained by numerical integration

        // need to get the reorganization energy
        double reorganization = 0.0;
        storage<double>* sdset = spectral_density->DirectAccess();
        for(int ind = 1; ind<numv; ind++)
            reorganization += sdset->data1D[ind]/ind;
        reorganization /= constants::pi;

        // doing numerical integrations (only for positve time)
        // putting initial values
        funa1[0]=0.0;
        funa2[0]= reorganization;
        for(int ind = 1; ind<numv; ind++)
        {
            double idel = 1.0/(deltat*deltat);
            complexv c0 = funa3[ind+1]*0.5*idel*times[ind-1]*times[ind];
            c0 += -funa3[ind]*idel*times[ind-1]*times[ind+1];
            c0 += funa3[ind-1]*0.5*idel*times[ind]*times[ind+1];

            complexv c1 = -funa3[ind+1]*0.5*idel*(times[ind-1]+times[ind]);
            c1 += funa3[ind]*idel*(times[ind-1]+times[ind+1]);
            c1 += -funa3[ind-1]*0.5*idel*(times[ind]+times[ind+1]);

            complexv c2 = funa3[ind+1]*0.5*idel;
            c2 += -funa3[ind]*idel;
            c2 += funa3[ind-1]*0.5*idel;

            funa2[ind]=funa2[ind-1]-coni*(c0*deltat+c1*0.5*(times[ind]*times[ind]-times[ind-1]*times[ind-1])+c2/3.0*(times[ind]*times[ind]*times[ind]-times[ind-1]*times[ind-1]*times[ind-1]) );

        }


        /// Im a2(inf) cannot be positive. This could happen because of integration errors
        // This cannot be allowed.
        // removing this part:
        if(funa2[numv-1].imag()>0)
        {
        	complexv shift = funa2[numv-1];

            for(int ind = 1; ind<numv; ind++)
            {
                funa2[ind] -= shift;
            }

            cout<<"Warning: notice that there is an accumulated error of a2(t) function:\n";
            cout<<"Im a2(inf) cannot be positive:\n";
            cout<<"The calculated value is: "<<shift<<":\n";
            cout<<"Whole a2(t) function has been shifted.\n";
        }


        for(int ind = 1; ind<numv; ind++)
        {
            double idel = 1.0/(deltat*deltat);
            complexv c0 = funa2[ind+1]*0.5*idel*times[ind-1]*times[ind];
            c0 += -funa2[ind]*idel*times[ind-1]*times[ind+1];
            c0 += funa2[ind-1]*0.5*idel*times[ind]*times[ind+1];

            complexv c1 = -funa2[ind+1]*0.5*idel*(times[ind-1]+times[ind]);
            c1 += funa2[ind]*idel*(times[ind-1]+times[ind+1]);
            c1 += -funa2[ind-1]*0.5*idel*(times[ind]+times[ind+1]);

            complexv c2 = funa2[ind+1]*0.5*idel;
            c2 += -funa2[ind]*idel;
            c2 += funa2[ind-1]*0.5*idel;

            funa1[ind]=funa1[ind-1]-coni*(c0*deltat+c1*0.5*(times[ind]*times[ind]-times[ind-1]*times[ind-1])+c2/3.0*(times[ind]*times[ind]*times[ind]-times[ind-1]*times[ind-1]*times[ind-1]) );


        }

    }




	// temporary functions
	asymptoticLF<double>* rep;// = accessory_a1P.DirectAccessRE();
	//asymptoticLF<double>* imp;// = accessory_a1P.DirectAccessIM();

    storage<complexv> sfuna(1);
    sfuna.Allocate(numv);


    // saving the calculated function A1: only positive part numv points
    // storing positive half
    for(int id = 0; id<numv; id++)
        sfuna.data1D[id] = funa1[id];
    *accessory_a1P = asymptoticLF_complexv(sfuna,0,deltat);
    // additional settings
	// rep->a = -ksi*temperature/const5305; //(funa1[numv+1].real()-funa1[numv].real())/deltat;
	rep = accessory_a1P->DirectAccessRE();
	rep->a = (funa1[numv-1].real()-funa1[numv-2].real())/deltat;
	rep->b = funa1[numv-1].real()+rep->a*deltat;
    rep = accessory_a1P->DirectAccessIM();
	rep->a = (funa1[numv-1].imag()-funa1[numv-2].imag())/deltat;
	rep->b = funa1[numv-1].imag()+rep->a*deltat;
	delete[] funa1;
    accessory_a1P->SetCausality(2,3);
    //cout<<"acc: "<<numv<<" "<<deltat<<"\n";
    //cout<<"accf:"<<accessory_a1P->GetXF()<<"\n";

	// saving the calculated function A2
    for(int id = 0; id<numv; id++)
        sfuna.data1D[id] = funa2[id];
    *accessory_a2P = asymptoticLF_complexv(sfuna,0,deltat);
    // additional settings
	rep = accessory_a2P->DirectAccessRE();
	rep->a = (funa2[numv-1].real()-funa2[numv-2].real())/deltat;
	rep->b = funa2[numv-1].real()+rep->a*deltat;;
    rep = accessory_a2P->DirectAccessIM();
	rep->a = (funa2[numv-1].imag()-funa2[numv-2].imag())/deltat;
	rep->b = funa2[numv-1].imag()+rep->a*deltat;;
	delete[] funa2;
    accessory_a2P->SetCausality(2,3);


	// saving the calculated function A3: only positive part
	// for negative part it is symmetric conjugate. Using that to improve accuracy.
    funa3[0] = funa3[0].real();
    for(int id = 0; id<numv; id++)
        sfuna.data1D[id] = funa3[id];
    *accessory_a3P = asymptoticLF_complexv(sfuna,0, deltat);
	rep = accessory_a3P->DirectAccessRE();
	rep->a = 0.0;//(funa3[numv+1].real()-funa3[numv].real())/deltat;
	rep->b = 0.0;//funa3[numv].real();
    rep = accessory_a3P->DirectAccessIM();
	rep->a = 0.0;//(funa3[numv+1].imag()-funa3[numv].imag())/deltat;
	rep->b = 0.0;//funa3[numv].imag();
	delete[] funa3;
    accessory_a3P->SetCausality(2,3);
    accessory_a3P->SetDiracAmplitude(globalgamma);

    //cout<<"acc3 "<< accessory_a3P->GetValueApproach(deltat);
    sfuna.Delete();

	delete[] dummySD;// = new double[numt];
	delete[] dummyPP;// = new double[numt];
    delete[] times;



    GenerateMBARFunction();

}



void numericalSD::GenerateAccessoryArrays(
double* dummySD, double* dummyPP,
complexv* funa1, complexv* funa2, complexv* funa3,
int numt,double deltav)
{

    // creates positive and negative parts consistent with FFT

	// at zero point adds small imaginary shift to keep convergence
	//complexv epsilon = complexv(0,deltav/50);

	int numt2 = numt/2;

	complexv* dataIn = new complexv[numt];
	complexv* dataOu = new complexv[numt];
	toolsFFT fplan;
        constants cc;

    if(funa1 != NULL){
        ////////////////////////////////////////////////////
	// calculating  A1 function: two parts; removing zero frequency point

        double dindv;

        dataIn[0] = 0.0;
        //dataIn[0] = dummySD[0]/epsilon/epsilon;

	// removing zero frequency point
	for(int indv = 1; indv < numt2; indv ++)
	{
                dindv = deltav*indv;
                dataIn[indv] = complexv( dummySD[indv]/dindv/dindv, 0.0);
                //cout<<dataIn[indv]<<" ---\n";
	}
	for(int indv = numt2; indv < numt; indv ++)
	{
                dindv = deltav*(indv-numt);
		dataIn[indv] = complexv( dummySD[indv]/dindv/dindv, 0.0);
                //cout<<dataIn[indv]<<" ---\n";
        }

        fplan.executeP(dataOu,dataIn,numt);

	// saving the calculated first part to function A1
	for(int indt = 0; indt < numt; indt ++)
        {
		funa1[indt] = dataOu[indt]*deltav/cc.pi2;
                //cout<<funa1[indt]<<"\n";

        }

	// calculating second part of A1 function
	dataIn[0] = 0.0;

	// removing zero frequency point
	for(int indv = 1; indv < numt2; indv ++)
	{
		dindv = deltav*indv;
                dataIn[indv] = complexv(dummyPP[indv]/dindv/dindv, 0.0);
	}
	for(int indv = numt2; indv < numt; indv ++)
	{
		dindv = deltav*(indv-numt);
		dataIn[indv] = complexv(dummyPP[indv]/dindv/dindv, 0.0);
	}
        fplan.executeP(dataOu,dataIn,numt);

	// saving the rest of calculated function A1
	for(int indt = 0; indt < numt; indt ++)
		funa1[indt] += dataOu[indt]*deltav/cc.pi2;


    }// done with function a1
        //cout<<funa1[0]<<"\n";
        //cout<<funa1[1]<<"\n";
        //cout<<funa1[2]<<"\n";

	if(funa2 != NULL){
	///////////////////////////////////////////////////////////
	// Now calculating function A2: this is done in two parts
        double dindv;

	// first removing area around zero; we could add this area late when saving the integral
	dataIn[0] = 0.0;
	for(int indv = 1; indv < numt2; indv ++)
	{
		dindv = deltav*indv;
		dataIn[indv] = complexv(dummySD[indv]/dindv, 0.0);
	}
	for(int indv = numt2; indv < numt; indv ++)
	{
		dindv = deltav*(indv-numt);
		dataIn[indv] = complexv(dummySD[indv]/dindv, 0.0);
	}
        fplan.executeP(dataOu,dataIn,numt);

        // saving  data
	for(int indt = 0; indt < numt; indt ++)
		funa2[indt] = dataOu[indt]*deltav/cc.pi2;


        double ksi = dummyPP[1]/deltav;
	dataIn[0] = complexv(ksi,0.0);
	for(int indv = 1; indv < numt2; indv ++)
	{
        	dindv = deltav*indv;
		dataIn[indv] = complexv( dummyPP[indv]/dindv, 0.0);
	}
	for(int indv = numt2; indv < numt; indv ++)
	{
		dindv = deltav*(indv-numt);
		dataIn[indv] = complexv( dummyPP[indv]/dindv, 0.0);
	}
        fplan.executeP(dataOu,dataIn,numt);

        // saving the remaining calculated function A2
        for(int indt = 0; indt < numt; indt ++)
                funa2[indt] += dataOu[indt]*deltav/cc.pi2;

	}// done with function a2

	if(funa3 != NULL){
                /////////////////////////////////////////////////////////////////////
                // Now calculating function A3
                for(int indv = 0; indv < numt; indv ++)
                	dataIn[indv] = complexv(dummySD[indv], 0.0);//+dummyPP[indv];
                fplan.executeP(dataOu,dataIn,numt);

                for(int indt = 0; indt < numt; indt ++)
                        funa3[indt] = dataOu[indt]*deltav/cc.pi2;

                for(int indv = 0; indv < numt; indv ++)
                {
                	dataIn[indv] = complexv(dummyPP[indv], 0.0);
                }
                fplan.executeP(dataOu,dataIn,numt);

                for(int indt = 0; indt < numt; indt ++)
                        funa3[indt] += dataOu[indt]*deltav/cc.pi2;

                //cout<<funa3[0]<<"\n";
	}// done with function a3

	delete[] dataIn;
	delete[] dataOu;

}

void numericalSD::GenerateMBARFunction()
{

    // doing numerical integrations

    // frequency parameters
    int numv = 2*spectral_density->GetN();
    double vmax =  2*spectral_density->GetXF();//    deltav*numv;
    double deltav = vmax/numv;


    storage<double> td(1);
    td.Allocate(numv);

//    double dummy1 = 1.0/2.0/temperature;
//    for(int find = 0; find<numv; find++)
//    {
//        double omegac = -deltav*numv/2 + deltav*find;
//
//        // now this funny summation
//        int summation = 100;
//        double sumerror = 1e64;
//        double sumval = 0.0;
//        int indn = 1;
//        while (fabs(sumerror)>fabs(sumval/10000.0) || indn <1000 )
//        {
//            double oml = omegac-indn*deltav;
//            double omh = omegac+indn*deltav;
//            double cotl = 0.0;
//            double coth = 0.0;
//            if(fabs(oml)>constants::smallepsilon)
//                cotl = (1.0+1.0/tanh(oml*dummy1))*spectral_density->Get(oml);
//            else
//                cotl = 2.0*temperature/deltav*spectral_density->Get(deltav);
//            if(fabs(omh)>constants::smallepsilon)
//                coth = (1.0+1.0/tanh(omh*dummy1))*spectral_density->Get(omh);
//            else
//                coth = 2.0*temperature/deltav*spectral_density->Get(deltav);
//
//            double indninv = 1.0/(double)indn;
//            sumerror = indninv*(-coth+cotl);
//            sumval += sumerror;
//
//            indn ++;
//
//        }
//
//        td.data1D[find] = sumval;
//    }
//    M_BAR_im =  asymptoticLF<double>(td,-deltav*numv/2,deltav, 0, 0, 0);
//
//    return;
//
    // otherwise
    // this must be done using Fourier transformations

	// M function is just a positive time Fourier transform of the accessory function A3
    // since a3 is the time domain correlation function

	// using accessory functions if available (this is one of the most accurate and fastest numerical results


    // the FFT dataset must be twice the length of the correlation function
    // since only the positive half of the correlation function is stored

    // FFT
    // time parameters
	int numt = numv;
    double deltat = const3000/vmax; // fs  = 2*pi*5305
	double tmax = numt*deltat; //

	// creating arrays for Fourier transformation
	complexv* dataIn = new complexv[numt];
	complexv* dataOu = new complexv[numt];
	toolsFFT fplan;

        // creating the dataset for FFT

    // positive time part
	for(int indt = 0; indt < numt/2; indt ++)
	{
		dataIn[indt] = accessory_a3P->Get(indt*deltat);
	}

    // negative time part is zero
	for(int indt = numt/2; indt < numt; indt ++)
		dataIn[indt] = 0.0;

	// zero point division due to half-Fourier
	dataIn[0] = complexv(dataIn[0].real()/2.0, 0.0);

    fplan.executeN(dataOu,dataIn,numt);

	// saving result only the imaginary part
    // positive
	for(int indt = 0; indt < numv/2; indt ++)
		td.data1D[indt+numv/2] = dataOu[indt].imag()*deltat;
    // negative
    for(int indt = numv/2; indt < numv; indt ++)
        td.data1D[indt-numv/2] = dataOu[indt].imag()*deltat;

    M_BAR_im =  asymptoticLF<double>(td,-deltav*numv/2,deltav, 0, 0, 0);

	delete[] dataIn;
	delete[] dataOu;

}


void numericalSD::CreateSymmetricFDCorrelation(double* sfSSD,double* gfuSD,double deltav,int numt)
{
	// this function generates the symmetric part of the frequency domain
	// correlation function from numerical spectral density gfuSD

	double dummy1 = 1.0/2.0/temperature;
        int numt2 = numt/2;

	for(int indv = 1; indv < numt2; indv ++)
	{
		double fre = indv*deltav;
		sfSSD[indv] = (1.0/tanh(fre*dummy1))*gfuSD[indv];
	}
	for(int indv = numt2; indv < numt; indv ++)
	{
		double fre = (indv-numt)*deltav;
		sfSSD[indv] = (1.0/tanh(fre*dummy1))*gfuSD[indv];
	}
	// doing quantum correction for zero frequency
	sfSSD[0] = 2.0*temperature*gfuSD[1]/deltav;
}

void numericalSD::CreateFDCorrelationSD(double* sfSSD,double* gfuSD,double deltav,int numt)
{
	// this function creates the frequency domain
	// correlation function from some form of spectral density
	// sfSSD is the even part of the correlation function
	// and gfuSD is the odd part (numerical spectral density)

	double dummy1 = 1.0/2.0/temperature;
        int numt2 = numt/2;

	for(int indv = 1; indv < numt2; indv ++)
	{
		double fre = indv*deltav;
		gfuSD[indv] = GetSD(fre);
                sfSSD[indv] = (1.0/tanh(fre*dummy1))*gfuSD[indv];
	}
	for(int indv = numt2+1; indv < numt; indv ++)
	{
		double fre = (indv-numt)*deltav;
		gfuSD[indv] = GetSD(fre);
                sfSSD[indv] = (1.0/tanh(fre*dummy1))*gfuSD[indv];
	}
	// doing quantum correction for zero frequency
	sfSSD[0] = 2.0*temperature*gfuSD[1]/deltav;
	gfuSD[0] = 0.0;
    sfSSD[numt2] = 0.0    ;
    gfuSD[numt2] = 0.0;
}


void numericalSD::printMBAR(string fname)
{
    int numv = spectral_density->GetN();
    double xmax = spectral_density->GetXF(); // cm-1
    double deltav = xmax/numv; // cm-1

    numv *=2;

    ofstream* file_str = new ofstream(fname.c_str());

    if(!file_str->is_open())
        cout<<"Error in numericalSD::printACCESSORY(string fname): "<<fname<<" file cannot been openned.\n";
    else
    {
        cout<<"Writing M(t) function\n";

        complexv val_s = 0.0;

        // in the file
				file_str->precision(16);

        *(file_str) << "# M(t), function of omega\n";

        for(int indd=0;indd<numv;indd++)
        {
            double time;
            time= -xmax+indd*deltav;

            complexv val_s = GetMf(time);
            *(file_str) << time;
            *(file_str) << "\t";
            *(file_str) << val_s.real();
            *(file_str) << "\t";

            *(file_str) << val_s.imag();
            *(file_str) << "\n";
        }
        file_str->close();
    }
    delete file_str;

}

void numericalSD::printACCESSORY(string fname)
{
	//here I must write a1(t), a2(t) and a3(t) function to a file in the same type format as the spectral density
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();


    asymptoticLF<double>* rep;// = accessory_a1P.DirectAccessRE();
	rep = accessory_a1P->DirectAccessRE();
	int numt = rep->GetN();
        //cout<<numt<<"\n";
	double tmax = rep->GetXF();
	double deltat = tmax/numt;

	ofstream* file_str = new ofstream(fname.c_str());

	if(!file_str->is_open())
		cout<<"Error in numericalSD::printACCESSORY(string fname): "<<fname<<" file cannot been openned.\n";
	else
	{
		cout<<"Writing a1(t), a2(t), a3(t) functions; units are forced to be the standard ones\n";

		complexv val_s = 0.0;

		// in the file
		// the 1 number is the number of data points
		// the 2 number is the frequency step
		// the data then follows

		file_str->precision(16);

		*(file_str) << "# a1(t), a2(t), a3(t) functions of time; units are forced to be the standard ones\n";

		for(int indd=0;indd<numt;indd++)
		{
			double time;
			time= indd*deltat;

            *(file_str) << time;
            *(file_str) << "\t";

            val_s = accessory_a1P->Get(time);
			*(file_str)<< val_s.real();
			*(file_str) << "\t";
			*(file_str)<< val_s.imag();
			*(file_str) << "\t";

			val_s = accessory_a2P->Get(time);
			*(file_str)<< val_s.real();
			*(file_str) << "\t";
			*(file_str)<< val_s.imag();
			*(file_str) << "\t";

			val_s = accessory_a3P->Get(time);
			*(file_str)<< val_s.real();
			*(file_str) << "\t";
			*(file_str)<< val_s.imag();
			*(file_str) << "\n";
		}
                file_str->close();
	}
        delete file_str;

}

void numericalSD::printGFUN(string fname)
{
	//here I must write g(t) function to a file in the same type format as the spectral density
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();

	if(lineshape == NULL)
		CreateLineshape();

    asymptoticLF<double>* rep;// = accessory_a1P.DirectAccessRE();
	rep = lineshape->DirectAccessRE();
	int numt = rep->GetN();
	double tmax = rep->GetXF();
	double deltat = tmax/numt;

        //cout<<numt<<"\n";
        //cout<<tmax<<"\n";
        //cout<<deltat<<"\n";


	ofstream* file_str = new ofstream(fname.c_str());

	if(!file_str->is_open())
		cout<<"Error in numericalSD::printGFUN(string fname): "<<fname<<" file cannot been openned.\n";
	else
	{
		cout<<"Writing g(t) function of time\n";

		complexv val_s = 0.0;

		// in the file
		// the 1 number is the number of data points
		// the 2 number is the frequency step
		// the data then follows

		file_str->precision(16);

		*(file_str) << "# g(t) function of time\n";

		for(int indd=0;indd<numt;indd++)
		{
			double time;
			time= indd*deltat;

            *(file_str)<< time*const5305;
            *(file_str) << "\t";

            val_s = GetGf(time);
			*(file_str)<< val_s.real();
			*(file_str) << "\t";
			*(file_str)<< val_s.imag();
			*(file_str) << "\n";
		}
                file_str->close();
	}
        delete file_str;

}

void numericalSD::printGDERIV1(string fname)
{
	//here I must write g'(t) function to a file in the same type format as the spectral density
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();


        asymptoticLF<double>* rep;// = accessory_a1P.DirectAccessRE();
	rep = accessory_a1P->DirectAccessRE();
	int numt = rep->GetN();
	double tmax = rep->GetXF();
	double deltat = tmax/numt;


	ofstream* file_str = new ofstream(fname.c_str());

	if(!file_str->is_open())
		cout<<"Error in numericalSD::printGDERIV1(string fname): "<<fname<<" file cannot been openned.\n";
	else
	{
		cout<<"Writing g'(t) functions\n";

		complexv val_s = 0.0;

		// in the file
		// the 1 number is the number of data points
		// the 2 number is the frequency step
		// the data then follows


		file_str->precision(16);
		*(file_str) << "# g'(t) function of time\n";

		for(int indd=0;indd<numt;indd++)
		{
			double time;
			time= indd*deltat;

            *(file_str)<< time*const5305;
            *(file_str) << "\t";

			val_s = GetGfD(time,1);
			*(file_str)<< val_s.real();
			*(file_str) << "\t";
			*(file_str)<< val_s.imag();
			*(file_str) << "\n";
		}
                file_str->close();
	}
        delete file_str;

}

void numericalSD::printGDERIV2(string fname)
{
	//here I must write g''(t) function to a file in the same type format as the spectral density

	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();

        asymptoticLF<double>* rep;// = accessory_a1P.DirectAccessRE();
	rep = accessory_a1P->DirectAccessRE();
	int numt = rep->GetN();
	double tmax = rep->GetXF();
	double deltat = tmax/numt;
        //cout<<deltat<<" test\n";



	ofstream* file_str = new ofstream(fname.c_str());

	if(!file_str->is_open())
		cout<<"Error in numericalSD::printGDERIV2(string fname): "<<fname<<" file cannot been openned.\n";
	else
	{
		cout<<"Writing g''(t)===C(t) function of time\n";

		complexv val_s = 0.0;

		// in the file
		// the 1 number is the number of data points
		// the 2 number is the frequency step
		// the data then follows


		file_str->precision(16);
		*(file_str) << "# g''(t)===C(t) function of time\n";

		for(int indd=0;indd<numt;indd++)
		{
			double time;
			time= indd*deltat;

            *(file_str)<< time*const5305;
            *(file_str) << "\t";

			val_s = accessory_a3P->Get(time);
			*(file_str)<< val_s.real();
			*(file_str) << "\t";
			*(file_str)<< val_s.imag();
			*(file_str) << "\n";
		}
                file_str->close();
	}
        delete file_str;

}

void numericalSD::printCORELF(string fname)
{
	printGDERIV2(fname);
}

void numericalSD::resolution(int& nump, double& freqstep)
{
	nump = spectral_density->GetN();
	freqstep = spectral_density->GetXF()/nump; // cm-1
}


// accessory functions
complexv numericalSD::GetA1(double tim)
{
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();
	return accessory_a1P->Get(tim/const5305);
}
complexv numericalSD::GetA2(double tim)
{
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();
	return accessory_a2P->Get(tim/const5305);
}
complexv numericalSD::GetA3(double tim)
{
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();
	return accessory_a3P->Get(tim/const5305);
}

complexv numericalSD::GetCf(double tim)
{
	return GetA3(tim);
}





asymptoticLF_complexv numericalSD::GetGf()
{
//	if(accessory_a1P == NULL)
//		ConstructAccessoryFunction();
	if(lineshape == NULL)
		CreateLineshape();
         //cout<<lineshape<<" test\n";

    asymptoticLF_complexv ret = *lineshape;
    ret.UpdateAxis(0.0,lineshape->GetStep()*const5305);

    return ret;
}


asymptoticLF_complexv numericalSD::GetGfD1f()
{
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();



//  coni*accessory_a2P->Get(tempt)
//			-coni*accessory_a2P->Get(0);

    asymptoticLF_complexv ret = (*accessory_a2P - accessory_a2P->Get(0) )*coni;

    ret.UpdateAxis(0.0,accessory_a3P->GetStep()*const5305);
	return ret;
}

asymptoticLF<double> numericalSD::GetSD()
{
//	if(accessory_a1P == NULL)
//		ConstructAccessoryFunction();
	if(spectral_density == NULL)
        {
		cout<<"ERROR: cannot return the spectral density\n";
         //cout<<lineshape<<" test\n";
                asymptoticLF<double> x;
                x = 0.0;
                return x;
        }

        return *spectral_density;
}

asymptoticLF_complexv numericalSD::GetCf()
{
	if(accessory_a1P == NULL)
		ConstructAccessoryFunction();

    asymptoticLF_complexv ret = *accessory_a3P;
    ret.UpdateAxis(0.0,accessory_a3P->GetStep()*const5305);
	return ret;
}

asymptoticLF_complexv numericalSD::GetMf()
{
    // this function returns the whole function
    // it returns in the same range as the spectral density itself.
    // however it makes positive and negative omegas so it creates
    // two times the number of points
    if(accessory_a1P == NULL)
        ConstructAccessoryFunction();



    int numv = spectral_density->GetN();
    double xmax = spectral_density->GetXF(); // cm-1
    double deltav = xmax/numv; // cm-1


    // making the dataset
    storage<complexv> idat(1);
    idat.Allocate(2*numv);

    // creating the dataset
    for(int iom = 0; iom<2*numv; iom ++)
    {
        double omega = -xmax +iom*deltav;

        idat.data1D[iom] = GetMf(omega);
    }

    //idat.data1D[0]=0.0;
    return asymptoticLF_complexv(idat,-xmax,deltav, 0.0,0.0,0.0);
}




asymptoticLF_complexv numericalSD::GetMfn()
{
    // this function returns the whole function
    // it returns in the same range as the spectral density itself.
    // however it makes positive and negative omegas so it creates
    // two times the number of points
    if(accessory_a1P == NULL)
        ConstructAccessoryFunction();



	int numv = spectral_density->GetN();
	double xmax = spectral_density->GetXF(); // cm-1
	double deltav = xmax/numv; // cm-1


    // making the dataset
    storage<complexv> idat(1);
    idat.Allocate(2*numv);

    // creating the dataset
    for(int iom = 0; iom<2*numv; iom ++)
    {
        double omega = -xmax +iom*deltav;

        idat.data1D[iom] = GetMfn(omega);
    }

    //idat.data1D[0]=0.0;
    return asymptoticLF_complexv(idat,-xmax,deltav, 0.0,0.0,0.0);
}


void numericalSD::SetT( double dtemperature)
{
    constants cc;

    temperature = const0_695*dtemperature;

}

void numericalSD::SetZPLGamma(double ig)
{
	globalgamma = ig;
}

void numericalSD::UpdateFlagAcc(int iflag)
{
	flagaccessories = iflag;
}
