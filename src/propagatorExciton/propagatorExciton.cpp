
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
    flagModred = 0; // 
    
    dmatrixT.SetDimension(2);
    dmatrix0.SetDimension(2);
    
    calcR = 0;
	superoperatorReady = 0;

	block=-1;
	numL=0;
	numR=0;
	shL=0;
	shR=0;

}

propagatorExciton::propagatorExciton(communicator_3rd& comm)
:propagatorMemory()
{
    meanEg = 0.0;
    meanEe = 0.0;
    meanEf = 0.0;

    flagMarkovian = 1;
    flagModred = 0; // 
	superoperatorReady = 0;

    dmatrixT.SetDimension(2);
    dmatrix0.SetDimension(2);
    
	block=-1;
	numL=0;
	numR=0;
	shL=0;
	shR=0;

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
	calcR->nonsecular = comm.nonsecular;

//	if(comm.nonsecular)
//        {
//	        calcR->AddMijkl(comm.zijkl);
//        }
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

	superoperatorReady = 1;
	if(flagMarkovian)
	{
		if(calcR->nonsecular)
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
	}
        else // this is always nonsecular
	{
		int numT=0;
		double deltaT=0;
        	superoperatorR = calcR->GetMemoryKernel(0, deltaT, numT , block);
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
        manifold0End = timesinternal.data1D[timeN-1]-constants::smallepsilon;

		//Initialize(superoperatorR, omegas_memory, dmatrixT, timesinternal)
	}
    }
    
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
		if(calcR->nonsecular){


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




