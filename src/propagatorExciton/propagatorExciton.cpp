
#include "propagatorExciton.hpp"
#include"../constructor-f-exciton/constructor-f-exciton.hpp"
#include"../constructor-f-exciton-gij/constructor-f-exciton-gij.hpp"
#include"../constants/constants.hpp"

propagatorExciton::propagatorExciton(int ig,int ie,int jf)
:propagatorMemory()
{
    numG = ig;
    numE = ie;
    numF = jf;

    meanEg = 0.0;
    meanEe = 0.0;
    meanEf = 0.0;


    flagMarkovian = 1; // Markovian
    flagModred = 0; // Population rates from modified Redfield
	flagNonsecular = 0; // secular only when markovian

    dmatrixT.SetDimension(2);
    dmatrix0.SetDimension(2);

    calcR = 0;
	superoperatorReady = 0;

	block=-1;
	numL=0;
	numR=0;
	shL=0;
	shR=0;

	// for multimode non-markovian relaxation
	DMMcontinous.SetDimension(1);

	DMMcontinousInternalTimeStep = 1.0;
    DMMcontinousInternalTimeNump = 0.0;

}

propagatorExciton::propagatorExciton(communicator_3rd& comm)
:propagatorMemory()
{
    meanEg = 0.0;
    meanEe = 0.0;
    meanEf = 0.0;

    flagMarkovian = 1;
    flagModred = 0; //
	flagNonsecular = 0; // secular only when markovian
	superoperatorReady = 0;

    dmatrixT.SetDimension(2);
    dmatrix0.SetDimension(2);

	block=-1;
	numL=0;
	numR=0;
	shL=0;
	shR=0;

	DMMcontinous.SetDimension(1);


    // the comm object must contain all excitonic information
    if(comm.excitonictransformed)
    {
        numG = comm.numG;
        numE = comm.numE;
        numF = comm.numF;


        SetEnergies(comm.evals);



    // specific part to the modified Redfield
    storage<asymptoticLF_complexv> g1un(1); // gfunctions first derivatives
    if(true){
    int numBosc = comm.mfun.GetSize();
    g1un.Allocate(numBosc);
    for(int ind=0; ind<numBosc; ind++)
    {
        numericalSD spd(comm.tempr,comm.spdf.data1D[ind]);
        g1un.data1D[ind] = spd.GetGfD1f();
    }
    }

    calcR = new calculator_redfield(energies, numG, numE, numF);
    calcR->AddFluctuations(comm.mfun,comm.cfun,comm.gfun,g1un,comm.kij,comm.gij,comm.zijkl);

//        calcR->AddFluctuations(comm.mfun,comm.cfun,comm.kij,comm.gij);
    calcR->tempr=comm.tempr;
	//calcR->nonsecular = comm.nonsecular;
	//flagNonsecular = comm.nonsecular;

    DMMcontinousInternalTimeStep= comm.propagationparameterStep;
    DMMcontinousInternalTimeNump= comm.propagationparameterNump;

    }
    else
    {
        std::cout<<"Error: propagatorExciton::propagatorExciton: transformation to eigenstates is necessary\n";
		calcR = 0;
    }

}


propagatorExciton::~propagatorExciton()
{
    if(calcR)
    {
        delete calcR;
        calcR = 0;
    }
	superoperatorReady = 0;

}


void propagatorExciton::SetEnergies(storage<double>& ien)
{
    energies = ien;


    meanEg = 0.0;
    meanEe = 0.0;
    meanEf = 0.0;

    for(int ind=0;ind<numG;ind++)
    {
        meanEg += ien.data1D[ind];
    }
    meanEg = meanEg /numG;
    for(int ind=0;ind<numE;ind++)
    {
        meanEe += ien.data1D[ind+numG];
    }
    meanEe = meanEe /numE;
    for(int ind=0;ind<numF;ind++)
    {
        meanEf += ien.data1D[ind+numG+numE];
    }
    meanEf = meanEf /numF;

}


void propagatorExciton::SetFluctuationMatrices
(storage<double>& matr00,storage<double>& matr10,storage<double>& matr20,
        storage<double>& matr11,storage<double>& matr21,storage<double>& matr22)
{
//    if(fAmplitudes.IsAlloc()) // complete matrix
//            return;
//
//    if(!mfun.IsAlloc())
//    {
//        std::cout<<"Error: propagatorExciton::SetFluctuationMatrices: set fluctuation functions\n";
//        return;
//    }
//
//    fAmplitudes.SetDimension(5);
//    int numBO = mfun.GetSize();
//    fAmplitudes.Allocate(numG+numE+numF,numG+numE+numF,numG+numE+numF,numG+numE+numF,numBO);
//
//    for(int inz=0; inz<numBO; inz++)
//    {
//        // setting gg block:
//        for(int ina=0; ina<numG; ina++)
//        for(int inb=0; inb<numG; inb++)
//        for(int inc=0; inc<numG; inc++)
//        for(int ind=0; ind<numG; ind++)
//        {
//            fAmplitudes.data5D[ina][inb][inc][ind][inz] = matr00.data5D[inz][ina][inb][inc][ind];
//        }
//
//        // setting eg block:
//        for(int ina=0; ina<numE; ina++)
//        for(int inb=0; inb<numE; inb++)
//        for(int inc=0; inc<numG; inc++)
//        for(int ind=0; ind<numG; ind++)
//      {
//            fAmplitudes.data5D[ina+numG][inb+numG][inc][ind][inz] = matr10.data5D[inz][ina][inb][inc][ind];
//            fAmplitudes.data5D[inc][ind][ina+numG][inb+numG][inz] = matr10.data5D[inz][ina][inb][inc][ind];
//        }
//
//        // setting fg block:
//        for(int ina=0; ina<numF; ina++)
//        for(int inb=0; inb<numF; inb++)
//        for(int inc=0; inc<numG; inc++)
//        for(int ind=0; ind<numG; ind++)
//        {
//            fAmplitudes.data5D[ina+numG+numE][inb+numG+numE][inc][ind][inz] = matr20.data5D[inz][ina][inb][inc][ind];
//            fAmplitudes.data5D[inc][ind][ina+numG+numE][inb+numG+numE][inz] = matr20.data5D[inz][ina][inb][inc][ind];
//        }
//
//        // setting ee block:
//        for(int ina=0; ina<numE; ina++)
//        for(int inb=0; inb<numE; inb++)
//        for(int inc=0; inc<numE; inc++)
//        for(int ind=0; ind<numE; ind++)
//        {
//            fAmplitudes.data5D[ina+numG][inb+numG][inc+numG][ind+numG][inz] = matr11.data5D[inz][ina][inb][inc][ind];
//        }
//
//        // setting fe block:
//        for(int ina=0; ina<numF; ina++)
//        for(int inb=0; inb<numF; inb++)
//        for(int inc=0; inc<numE; inc++)
//        for(int ind=0; ind<numE; ind++)
//        {
//            fAmplitudes.data5D[ina+numG+numE][inb+numG+numE][inc+numG][ind+numG][inz] = matr21.data5D[inz][ina][inb][inc][ind];
//            fAmplitudes.data5D[inc+numG][ind+numG][ina+numG+numE][inb+numG+numE][inz] = matr21.data5D[inz][ina][inb][inc][ind];
//        }
//
//        // setting ff block:
//        for(int ina=0; ina<numF; ina++)
//        for(int inb=0; inb<numF; inb++)
//        for(int inc=0; inc<numF; inc++)
//        for(int ind=0; ind<numF; ind++)
//        {
//            fAmplitudes.data5D[ina+numG+numE][inb+numG+numE][inc+numG+numE][ind+numG+numE][inz] =
//                    matr22.data5D[inz][ina][inb][inc][ind];
//        }
//    }
}
void propagatorExciton::SetFluctuationFunctions
(storage<asymptoticLF_complexv>& imfun,storage<asymptoticLF_complexv>& icfun,double ttt)
{
//    mfun = imfun;
//    cfun = icfun;
//    tempr = ttt;
    if(calcR == 0)
    {
        calcR = new calculator_redfield(energies, numG, numE, numF);
    }
    calcR->mfun = imfun;
    calcR->cfun = icfun;
    calcR->tempr = ttt;
}

// sets Boltzmann equilibrium distribution
void propagatorExciton::MakeEquilibriumSim()
{
}

void propagatorExciton::SetCurrentDM(storage<complexv>& imtr)
{
    dmatrix0 = imtr;
}
void propagatorExciton::SetBlock(int iblock)
{
	block=iblock;
	if(block == 0)
	{
		numL=numG;
    		numR=numG;
    		shL = 0;
    		shR = 0;
	}
	else if(block == 10)
	{
		numL=numE;
    		numR=numG;
    		shL = numG;
    		shR = 0;
	}
	else if(block == 11)
	{
		numL=numE;
    		numR=numE;
    		shL = numG;
    		shR = numG;
	}
	else if(block == 20)
	{
		numL=numF;
    		numR=numG;
    		shL = numG+numE;
    		shR = 0;
	}
	else if(block == 21)
	{
		numL=numF;
    		numR=numE;
    		shL = numG+numE;
    		shR = numG;
	}
	else if(block == 22)
	{
		numL=numF;
    		numR=numF;
    		shL = numG+numE;
    		shR = numG+numE;
	}
	else
	{
		std::cout<<"Error: wrong block label in propagatorExciton::SetBlock(int block)\n";
	}
}

storage<complexv> propagatorExciton::PropagateDM(storage<double>& times)
{
    	if (!superoperatorReady)
    	{
				if(calcR==0)
				{
					cout<<"Error: calculator in propagatorExciton::PropagateDMB is not ready\n";
				}
        energiesReorg =  calcR->GetReorganizations();

				if(flagMarkovian)
				{
					if( flagNonsecular )
	            superoperatorM = calcR->GetRelaxationSuperoperatorM(block);
					else
					{
            if(flagModred)
                superoperatorSP = calcR->GetTransportRatesModRed();
            else
                superoperatorSP = calcR->AddTransportRates();
						calcR->AddLifetimeDephasings();
						superoperatorSG = calcR->AddPureDephasings();
					}
					superoperatorReady = 1;
				}
				// now nonMarkovian
        else if(flagMemoryWithCfun) // this is always nonsecular
				{ // this is  memory relaxation kernel
					// factorized into amplitudes (kernel)
					// and correlation functions
					cout<<"Preparing relaxation amplitudes for non-Markovian propagation\n";

					// first about times:
					// int numT=0;
		      // double deltaT=0;
		      // if (internaltimeN != 0)
		      // {
		      //   numT=internaltimeN;
		  		// 	deltaT=internaltimeS;
		      // }
					// else
					// {
					// 	// define history parameters from the zeroth correlation function
					// 	deltaT = calcR->cfun.data1D[0].GetStep();
					// 	numT = calcR->cfun.data1D[0].GetN();
					// }
					// DMMcontinousInternalTimeStep = deltaT;
			    // DMMcontinousInternalTimeNump = numT;
					// cout << "internal time points: "<< numT << "\n";
					// cout << "internal time steps: " << deltaT << "\n";


					superoperatorR = calcR->GetMemoryKern(block);
	        kernel = &superoperatorR;
        	DMM = &dmatrixT;

					// setting up omegas:
					omegas_memory.Allocate(numL*numR);
					for(int id2=0;id2<numL;id2++)
					for(int id1=0;id1<numR;id1++)
					{
						omegas_memory.data1D[id2*numR+id1]=energies.data1D[id2+shL]-energies.data1D[id1+shR];
					}
					// also setting the history of density matrix into 1D form:
					//double DMMcontinousInternalTimeStep;
	        //int    DMMcontinousInternalTimeNump;
					// must be defined
					if(!dmatrixT.IsAlloc())
					{
						dmatrixT.Allocate(DMMcontinousInternalTimeNump,numL*numR);
						for(int ind2=0;ind2<numL;ind2++)
						for(int ind1=0;ind1<numR;ind1++)
							dmatrixT.data2D[0][ind2*numR+ind1]=dmatrix0.data2D[ind2][ind1];
						// filling up times array
						timesinternal.Delete();
						timesinternal.Allocate(DMMcontinousInternalTimeNump);
						timesinternal.FillLinear(0,-DMMcontinousInternalTimeStep,DMMcontinousInternalTimeNump);
					}
					superoperatorReady = 1;

				}
				else
				{
					// this is full memory relaxation kernel as a function of time
					cout<<"Preparing complete non-Markovian relaxation tensor\n";
      		int numT=0;
      		double deltaT=0;
      		if (internaltimeN != 0)
      		{
        		numT=internaltimeN;
  					deltaT=internaltimeS;
        		cout << "internal time points: "<< numT << "\n";
        		cout << "internal time steps: " << deltaT << "\n";
      		}
      		superoperatorR = calcR->GetMemoryKernel(0, deltaT, numT , block);
					DMMcontinousInternalTimeStep = deltaT;
	    		DMMcontinousInternalTimeNump = numT;
					//cout<<"DMMcontinousInternalTimeStep = "<<DMMcontinousInternalTimeStep<<"\n";
					//cout<<"DMMcontinousInternalTimeNump = "<<DMMcontinousInternalTimeNump<<"\n";
					// converting into lower dimensional form
					superoperatorR = Convert5DTo3D(superoperatorR);
					// getting sizes
					int num2,num1;
					superoperatorR.GetSize(numT,num2,num1);

					// setting up omegas:
					omegas_memory.Allocate(numL*numR);
					for(int id2=0;id2<numL;id2++)
					for(int id1=0;id1<numR;id1++)
					{
						omegas_memory.data1D[id2*numR+id1]=energies.data1D[id2+shL]-energies.data1D[id1+shR];
					}

					// also setting the history of density matrix into 1D form:
					if(!dmatrixT.IsAlloc())
					{
						dmatrixT.Allocate(numT,numL*numR);
						for(int ind2=0;ind2<numL;ind2++)
						for(int ind1=0;ind1<numR;ind1++)
							dmatrixT.data2D[0][ind2*numR+ind1]=dmatrix0.data2D[ind2][ind1];

						// filling up times array
						timesinternal.Delete();
						timesinternal.Allocate(numT);
						timesinternal.FillLinear(0,-deltaT,numT);
					}

        		kernel = &superoperatorR;
        		DMM = &dmatrixT;
        		int timeN = timesinternal.GetSize();
        	//manifold0End = timesinternal.data1D[timeN-1]-constants::smallepsilon;
					superoperatorReady = 1;
				}

				cout << "Relaxation kernel prepared.\n";
				//exit(0);

    }// done with preparing the kernels

    storage<complexv> result(3);

    // getting timestep
    if(!times.IsAlloc())
    {
        cout<<"Error: input argument must be a time series (at leat two elements)\n";
        return result;
    }
    if(times.GetSize()==1)
    {
        cout<<"Error: input argument must be a time series (at leat two elements)\n";
        return result;
    }

    double td = times.data1D[1]-times.data1D[0];
    double ti = times.data1D[0];
    int tvals = times.GetSize();
    result.Allocate(tvals,numL,numR);

    // propagation from time0 to timet;
    if(flagMarkovian)
    {
		if(flagNonsecular){
			cout<<"Propagating Markovian NonSecular\n";

	if(pevals == 0)
	{
		// making Liouville superoperator
	        storage<complexv> matrix(2);
	        matrix.Allocate(numL*numR,numL*numR);

	        for(int il2 = 0; il2<numL;il2++)
        	for(int il1 = 0; il1<numR;il1++)
        	{
        	    int il = il2*numR+il1;
	            for(int ir2 = 0; ir2<numL;ir2++)
	            for(int ir1 = 0; ir1<numR;ir1++)
	            {
	                int ir = ir2*numR+ir1;

	                matrix.data2D[il][ir] = 0.0;
	                if(il2 == ir2 && il1 == ir1){
	                    matrix.data2D[il][ir] +=  cnni*(energies.data1D[il2+shL]-energiesReorg.data2D[il2+shL][il2+shL]);
	                    matrix.data2D[il][ir] +=  coni*(energies.data1D[il1+shR]-energiesReorg.data2D[il1+shR][il1+shR]);
	                }
	                matrix.data2D[il][ir] += superoperatorM.data4D[il2][il1][ir2][ir1];
	            }
	        }
	        InitMaster(matrix);
	}

        // Markovian propagation
        // propagation is realized by using the diagonalization of relaxation superoperator
        storage<complexv> y = Convert2DTo1D(dmatrix0);
        storage<complexv> retval;

        for(int numt = 0; numt< tvals; numt++)
        {
            retval =  Get(y,times.data1D[numt]);

            for(int il2 = 0; il2<numL;il2++)
            for(int il1 = 0; il1<numR;il1++)
            {
                int il = il2*numR+il1;
                result.data3D[numt][il2][il1] = retval.data1D[il];
                if(numt == tvals-1)
                    dmatrix0.data2D[il2][il1] = retval.data1D[il];
            }
        }
		}
		else
		{
			// secular dynamics

			cout<<"Propagating Markovian Secular\n";
	// part A - coherences
        // propagation is realized by using direct exponentiation
	        for(int numt = 0; numt< tvals; numt++)
	        {
			double time;
            		for(int il2 = 0; il2<numL;il2++)
        		for(int il1 = 0; il1<numR;il1++)
			if((il2+shL)!=(il1+shR))
            		{
                          complexv argument = 0;
                          argument += cnni*(energies.data1D[il2+shL]-energiesReorg.data2D[il2+shL][il2+shL]);
                          argument += coni*(energies.data1D[il1+shR]-energiesReorg.data2D[il1+shR][il1+shR]);
                          argument -= superoperatorSG.data2D[il2+shL][il1+shR];
	                  result.data3D[numt][il2][il1] = exp(argument*times.data1D[numt])*dmatrix0.data2D[il2][il1];

	                  if(numt == tvals-1)
        	            dmatrix0.data2D[il2][il1] = result.data3D[numt][il2][il1];
		        }
        	}

	if(block ==0 ||block ==11 ||block ==22)
	{
		// part B - populations

		if(pevals == 0)
		{
	        	storage<double> matrix(2);
		        matrix.Allocate(numL,numL);
	        	for(int il2 = 0; il2<numL;il2++)
			for(int ir2 = 0; ir2<numL;ir2++)
	            	{
	                	matrix.data2D[il2][ir2] = superoperatorSP.data2D[il2+shL][ir2+shL];
	            	}
		        InitMaster(matrix);
		}
        	storage<double> y(1);
        	y.Allocate(numL);
        	for(int il2 = 0; il2<numL;il2++)
            		y.data1D[il2] = dmatrix0.data2D[il2][il2].real();
	        storage<double> retval;

	        // propagation is realized by using the diagonalization of relaxation superoperator
	        for(int numt = 0; numt< tvals; numt++)
	        {
	            retval =  Get(y,times.data1D[numt]);

            		for(int il2 = 0; il2<numL;il2++)
            		{
	                  result.data3D[numt][il2][il2] = retval.data1D[il2];

	                  if(numt == tvals-1)
        	            dmatrix0.data2D[il2][il2] = retval.data1D[il2];

		        }
        	}
	}}
    }
    else
    {
	// this is non-markovian part
        //cout<<"Error: sorry, this part of code is missing\n";
	// I will propagate the Interaction representation internally

	cout<<"Propagating NonMarkovian\n";
	result =  Propagate(times);
	result = Convert2DTo3D(result,numL,numR);

            		for(int il2 = 0; il2<numL;il2++)
            		for(int il1 = 0; il1<numR;il1++)
        	            dmatrix0.data2D[il2][il1] = result.data3D[tvals-1][il2][il1];

    }
    return result;
}
void propagatorExciton::SetMarkovian(int iflag)
{
    flagMarkovian = iflag;
}


void propagatorExciton::SetDeltaDM(int& ie,int &ig)
{
	if(dmatrixT.IsAlloc())
	{
		int nl,nt;
		dmatrixT.GetSize(nt,nl);
		if(nl!= numL*numR)
		{
			dmatrixT.Delete();
			dmatrixT.Allocate(nt,numL*numR);
		}


    		for(int it = 0; it< nt; it++)
    		for(int ig1 = 0; ig1< numL*numR; ig1++)
    		{
        		dmatrixT.data2D[it][ig1] = 0.0;
    		}
    		dmatrixT.data2D[0][ie*numR+ig] = 1.0;

		// filling up times array
		double deltaT=timesinternal.data1D[0]-timesinternal.data1D[1];
		timesinternal.FillLinear(0,-deltaT,nt);

	}

    if(!dmatrix0.IsAlloc())
	dmatrix0.Allocate(numL,numR);
    else
	{
		int nl,nr;
		dmatrix0.GetSize(nl,nr);
		if(nl!= numL || nr!= numR)
		{
			dmatrix0.Delete();
			dmatrix0.Allocate(numL,numR);
		}
	}


    for(int ig1 = 0; ig1< numL; ig1++)
    for(int ig2 = 0; ig2< numR; ig2++)
    {
        dmatrix0.data2D[ig1][ig2] = 0.0;
    }
    dmatrix0.data2D[ie][ig] = 1.0;

}



// the following functions upgrate the propagatorMemory with
// explicit block form for the excitons

void propagatorExciton::ConvoluteGen(
    storage<complexv>& der  // current derivative values
)
{
    if(flagMemoryWithCfun)
        ConvoluteCfun(der);
    else
        Convolute(der);  // original approach in propagatorMemory
}

void propagatorExciton::ConvoluteCfun(
    storage<complexv>& der  // current derivative values
)
{
	// timesinternal are internal time parameters for the memory kernel
	// its timestep is equal to DMMcontinousInternalTimeStep
    double timeini = timesinternal.data1D[0];
    int itemp,dimension;
    DMM->GetSize(itemp,dimension);
    int numosc = calcR->cfun.GetSize();


    for(int ia=0;ia<numL;ia++)
    for(int ib=0;ib<numR;ib++)
    {
    	der.data1D[ia*numR+ib] = 0.0;
    	double omegalxtimel = (energies.data1D[ia+shL]-energies.data1D[ib+shR])*timeini;

 			for(int io=0;io<numosc;io++)
    	{
					// these time parameters are for integrating the kernel
					double deltat = calcR->cfun.data1D[io].GetStep();
					double numT = calcR->cfun.data1D[io].GetN();
					// unless specified externally:
		      if (internaltimeN != 0)
		      {
						numT=internaltimeN;
		  			deltat=internaltimeS;
		      }
					DMMcontinousInternalTimeStep = deltat;
			    DMMcontinousInternalTimeNump = numT;
					cout << "internal time points: "<< numT << "\n";
					cout << "internal time steps: " << deltat << "\n";


    			for(int indt=0;indt<numT; indt++)
    			{
	    			// notice: "time" goes in reverse direction
    				double time = timeini - deltat*indt;

        		// this "tau" goes forward direction starting from zero
        		double tau = indt*deltat;

						complexv cfun_requested = calcR->cfun.data1D[io].Get(tau);

						// summation over internal kernel indices
        		for(int ia1=0;ia1<numL;ia1++)
        		for(int ib1=0;ib1<numR;ib1++)
    				{
        			// recalculating the relaxation tensor - copy from calculator_redfield
							complexv R1=0;
							if(ib==ib1)
	    				for(int ic=0;ic<numL;ic++)
							{
								double omegat = energies.data1D[ic+shL]-energies.data1D[ib1+shR];
								R1 += exp(cnni*omegat*tau)*kernel->data3D[io][(ia+shL)*numL+(ic+shL)][(ic+shL)*numL+(ia1+shL)]
                    *cfun_requested;
							}
							complexv R4=0;
							if(ia==ia1)
        			for(int id=0;id<numR;id++)
							{
								double omegat = energies.data1D[ia+shL]-energies.data1D[id+shR];
								R4 += exp(cnni*omegat*tau)*conj(kernel->data3D[io][(ib1+shR)*numR+(id+shR)][(id+shR)*numR+(ib+shR)]
                    *cfun_requested);
							}
							complexv R2=0;
							{
								double omegat = energies.data1D[ia1+shL]-energies.data1D[ib+shR];
								R2 += exp(cnni*omegat*tau)*conj(kernel->data3D[io][(ia+shL)*numL+(ia1+shL)][(ib1+shR)*numR+(ib+shR)]
                    *cfun_requested);
							}
							complexv R3=0;
							{
								double omegat = energies.data1D[ia+shL]-energies.data1D[ib1+shR];
								R3 += exp(cnni*omegat*tau)*kernel->data3D[io][(ia+shL)*numL+(ia1+shL)][(ib1+shR)*numR+(ib+shR)]
                    *cfun_requested;
							}

							// now making the interaction picture transformation
							double omega = energies.data1D[ia1+shL]-energies.data1D[ib1+shR];

        			der.data1D[ia*numR+ib] =  -(R1-R2-R3+R4)*exp(-coni*(omegalxtimel-omega*time))
                *deltat* DMMcontinous.data1D[ia1*numR+ib1].Get(tau);
      		}// internal indices a1 b1
    		}// summation over time


			}// summation over bath oscillators

    }// summation over ia ib
} // end of function



// propagates DMM with respect to itself and history DMM1 and DMM2
storage<complexv> propagatorExciton::Propagate(storage<double>& times)
{

    // external propagation parameters
    double tdext = times.data1D[1]-times.data1D[0];
    double tiext = times.data1D[0];
    int tvalsext = times.GetSize();

		//cout<<"tdext "<<tdext<<"\n";
		//cout<<"tiext "<<tiext<<"\n";
		//cout<<"tvalsext "<<tvalsext<<"\n";

    // internal propagation parameters
    double tdint = timesinternal.data1D[0]-timesinternal.data1D[1];
    int numT,dimension;
    DMM->GetSize(numT,dimension);

	if(flagMemoryWithCfun)
	// create continous representation of DMM
	// with time starting from zero and going upwards
	{
		if(DMMcontinous.IsAlloc()) DMMcontinous.Delete();
		DMMcontinous.Allocate(dimension);
		for(int ix=0;ix<dimension;ix++)
		{
		storage<complexv> dataset(1);
		dataset.Allocate(numT+2);
		for(int it=0; it<numT; it++)
			dataset.data1D[it]=DMM->data2D[it][ix];
		dataset.data1D[numT]=0.0;
		dataset.data1D[numT+1]=0.0;
		asymptoticLF_complexv tfun(dataset,0.0,tdint);
		DMMcontinous.data1D[ix]=tfun;
		}
	}


    storage<complexv> ret(2);
    ret.Allocate(tvalsext,dimension);

    storage<complexv> derivs(1);
    derivs.Allocate(dimension);

    //loops over external times
    for(int itext=0;itext<tvalsext;itext++)
    {
        double timeext = times.data1D[itext];

        // propagation with small internal timesteps;
				double tiint = timesinternal.data1D[0];
        int numloops = (int)((timeext-tiint)/tdint); // this should return floor value
        double reminder = timeext-tiint-numloops*tdint;

				//cout<<"space\n";
				//cout<<"timeext = "<<timeext<<"\n";
				//cout<<"numloops = "<<numloops<<"\n";
				//cout<<"reminder = "<<reminder<<"\n";

        for(int it=0;it<numloops;it++)
        {
            ConvoluteGen( derivs);
            Update(derivs,tdint); // tdint is optional
        }


        // doing reminder assignment:
        ConvoluteGen( derivs);
        // and new result
				// DMM matrix is not updated
        for(int ap = 0; ap<dimension; ap++)
        {
            double omega = omegas_memory.data1D[ap];
            ret.data2D[itext][ap] = exp(cnni*omega*timeext)*(DMM->data2D[0][ap] + derivs.data1D[ap]*reminder);
        }
				//timesinternal.FillLinear(
				//		timesinternal.data1D[0]+reminder,
				//		timesinternal.data1D[1]-timesinternal.data1D[0],
				//		timesinternal.GetSize());

    }

    return ret;
}
