#include"toolsRandom.hpp"

#include"../complexv/complexv.hpp"
#include"../constants/constants.hpp"
#include"../toolsFFT/toolsFFT.hpp"

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

toolsRandom::toolsRandom()
{
    srand (time(NULL));
}

toolsRandom::toolsRandom(int ival)
{
    srand (ival);
}

double toolsRandom::RandLD(double k)//returns linear random from 0 to <k
{
  double s;

  s=(double)rand();
  //cout<<"\n"<<s%k<<"\n";
  return k*(s/RAND_MAX);
}


double toolsRandom::RandED(double k)//returns exponential random
{
  double s;

  s=RandLD(1);

  return -k*log(s);
}
double toolsRandom::RandGD(double v, double s)//returns gaussian random
{
  double x,y;

  do
    {
      y=RandLD(1);
      x=-4+RandLD(8);
    }
  while(y>exp(-x*x/2));
  return x*s+v;
}

int toolsRandom::GaussianTrajectory(complexv* gtraj, int nump, double correlationL)
{
	// this function returns the trajectory 
	// characterized by a Gaussian statistics and by an exponential correlation decay
	// this should be the Overdamped Brownian Oscillator model

	// correlationL is in the units of trajectory points
	// the width of Gaussian is unit.

	// the trajectory is created using the Fourier amplitude method

	complexv* tft = new complexv[nump];
        
	constants cst;
        const double const_2pi = cst.pi2;

	correlationL = 1.0/(correlationL/nump);

	for(int in=0;in<nump;in++)
	{
		if( in<nump/2 )
		{
			double kdk = const_2pi*in;
			tft[in] = nump * sqrt( 2.0*correlationL/(kdk*kdk + correlationL*correlationL) );

			// adding random phase
			double phase = RandLD(const_2pi);
			tft[in] *= complexv(cos(phase),sin(phase));
		}
		else
		{
			double knk = const_2pi*(in-nump);
			tft[in] = nump * sqrt( 2.0*correlationL/(knk*knk + correlationL*correlationL) );

			double phase = RandLD(const_2pi);
			tft[in] *= complexv(cos(phase),sin(phase));
		}
	}

	// now it is time to do FFT

	toolsFFT tfft;
        tfft.executeP(gtraj,tft,nump);

        return 0;
	// that is it ! DONE!
}



