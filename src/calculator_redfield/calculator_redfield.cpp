#include"../toolsFFT/toolsFFT.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../storage/storage.hpp"
#include"../constants/constants.hpp"
#include"../propagatorM/propagatorM.hpp"
#include"../numericalSD/numericalSD.hpp"

#include"calculator_redfield.hpp"
    
calculator_redfield::calculator_redfield()
{
    ready = 0;
    Setup();

}
calculator_redfield::calculator_redfield(
storage<double>& ielevs)
{

    Setup();

    // acquire system parameters
    evals = ielevs;
        
    ready = 0;
}

calculator_redfield::calculator_redfield(
storage<double>& ielevs, int numg, int nume, int numf)
{

    Setup();

    // acquire system parameters
    evals = ielevs;
        
    ready = 0;
    numG = numg;
    numE = nume;
    numF = numf;
}

void calculator_redfield::AddFluctuations
        (storage<asymptoticLF_complexv>& imfun,
         storage<asymptoticLF_complexv>& icfun,
         storage<asymptoticLF_complexv>& igfun,
         storage<asymptoticLF_complexv>& igdun,
         storage<double>& ikij,
         storage<double>& igij,
         storage<double>& iMijkl )
{
    cfun = icfun;
    gfun = igfun;
    gdun = igdun;
    Mijkl = iMijkl;
    AddFluctuations(imfun,ikij,igij);
}

void calculator_redfield::AddFluctuations
        (storage<asymptoticLF_complexv>& imfun,storage<double>& ikij,storage<double>& igij)
{
    mfun = imfun;
    kij = ikij; // amplitudes for all pairs
    gij = igij; // amplitudes for all pairs

    ready = 1;
    
    int numl = evals.GetSize();
    int numBosc = mfun.GetSize();


    // calculating reorganization energies
    // calculaTING REOGANIZATION functions
    if(!reorganizations.IsAlloc())
        reorganizations.Allocate(numl,numl);
    if(true)
    {
        cout<<"# reorganization energies:\n";

        // calculating the partition function
        for(int ind = 0; ind<numl; ind++)
            for(int inr = 0; inr<numl; inr++)
            {

                double part = 0;

                for(int inb = 0; inb<numBosc; inb++)
                {
                    part += gij.data3D[ind][inr][inb]*(mfun.data1D[inb].Get(0)).imag();
                }
                reorganizations.data2D[ind][inr] = -part;

                if(inr == numl-1)
                    cout<<reorganizations.data2D[ind][inr]<<"\n";
                else
                    cout<<reorganizations.data2D[ind][inr]<<"\t";
            }
    }
}

storage<double> calculator_redfield::GetReorganizations()
{
	if(!reorganizations.IsAlloc())
	{
		cout<<"Error: cannot return reorganizations\n";
	}
	return reorganizations;
}


storage<double> calculator_redfield::AddTransportRates()
{
    // here I create redfield rates and 
    // dephasing constants from mfun
    
    int numG = evals.GetSize();
    
    if(!rates.IsAlloc())
    {
        rates.Allocate(numG,numG);

        // number of bath oscillators
        int numBosc = mfun.GetSize();


        for(int iosc = 0; iosc<numBosc; iosc ++)
        {
        
        	for(int ia = 0; ia<numG; ia++)
        		for(int ib = 0; ib<numG; ib++)
        			if(ia != ib)
        			{
            	double omegaij = evals.data1D[ia]-evals.data1D[ib];
            	//omegaij += -reorganizations.data2D[ia][ia] + reorganizations.data2D[ib][ib];
                rates.data2D[ia][ib] += 2.0*kij.data3D[ia][ib][iosc]*(mfun.data1D[iosc].Get(-omegaij)).real();
        			}
        }
    
        // filling diagonal values
        for(int ia = 0; ia<numG; ia++)
        {
        	rates.data2D[ia][ia] = 0.0;
        	for(int ib = 0; ib<numG; ib++)
        		if(ia != ib )
        		{
        			rates.data2D[ia][ia] -= rates.data2D[ib][ia];
        		}
        }

        cout<<"# Redfield rate matrix:\n";
        for(int ia = 0; ia<numG; ia++)
        {
        	for(int ib = 0; ib<numG; ib++)
        	{
        		if(ib == numG-1)
        			cout<<rates.data2D[ia][ib]<<"\n";
        		else
        			cout<<rates.data2D[ia][ib]<<"\t";
        	}
        }
    }
    
    
    return rates;
    
}

storage<double> calculator_redfield::GetTransportRatesModRed()
{

    if(!Mijkl.IsAlloc())
    {
        cout<<"Error: calculator_redfield::GetRelaxationSuperoperatorM(): fluctuation matrix is missing\n";
        return rates;
    }

	complexv laccc = 0.0;
	complexv lcacc = 0.0;
	complexv lcccc = 0.0;
	complexv lccaa = 0.0;
	complexv Fcccc = 0.0;
	complexv Fccaa = 0.0;
	complexv Faaaa = 0.0;
	complexv Faacc = 0.0;
	complexv Fcaaa = 0.0;
	complexv Fcacc = 0.0;
	complexv Facaa = 0.0;
	complexv Faccc = 0.0;
	complexv Fcaac = 0.0;
    
    int num = evals.GetSize();
    
    if(!rates.IsAlloc())
    {
        rates.Allocate(num,num);

        // number of bath oscillators
        int numBosc = mfun.GetSize();

    

        	for(int ia = 0; ia<num; ia++)
        		for(int ic = 0; ic<num; ic++)
        			if(ia != ic)
        			{
                        
                        
    
                // updating reorganization energies from M 
                laccc = 0.0;
                lcacc = 0.0;
                lcccc = 0.0;
                lccaa = 0.0;
	
                for(int inb = 0; inb<numBosc; inb++)
                {
                    laccc += Mijkl.data5D[ia][ic][ic][ic][inb]*(mfun.data1D[inb].Get(0)).imag();
                    lcacc += Mijkl.data5D[ic][ia][ic][ic][inb]*(mfun.data1D[inb].Get(0)).imag();
                    lcccc += Mijkl.data5D[ic][ic][ic][ic][inb]*(mfun.data1D[inb].Get(0)).imag();
                    lccaa += Mijkl.data5D[ic][ic][ia][ia][inb]*(mfun.data1D[inb].Get(0)).imag();
                }
                laccc = -laccc;
                lcacc = -lcacc;
                lcccc = -lcccc;
                lccaa = -lccaa;

                
        complexv fval = 0.0;

        for(int inb = 0; inb<numBosc; inb ++)
        {
    
        // getting time integration parameters
        int numt = gfun.data1D[inb].GetN();
        double tmax = gfun.data1D[inb].GetXF();
        double step = tmax/numt; 


        // time integration
		for( int indin=0; indin < numt; indin++ )
		{
			double ti = indin*step;

			//getting required lineshape functions
			complexv g0 = gfun.data1D[inb].Get(ti);
			complexv g1 = gdun.data1D[inb].Get(ti);
			complexv g2 = cfun.data1D[inb].Get(ti);

			// g functions
			complexv G0Scccc = Mijkl.data5D[ic][ic][ic][ic][inb]*g0;
			complexv G0Sccaa = Mijkl.data5D[ic][ic][ia][ia][inb]*g0; 
			complexv G0Saaaa = Mijkl.data5D[ia][ia][ia][ia][inb]*g0; 
			complexv G0Saacc = Mijkl.data5D[ia][ia][ic][ic][inb]*g0; 

			// g 1 derivative functions
			complexv G1Saaac = Mijkl.data5D[ia][ia][ia][ic][inb]*g1; 
			complexv G1Sccac = Mijkl.data5D[ic][ic][ia][ic][inb]*g1; 
			complexv G1Scacc = Mijkl.data5D[ic][ia][ic][ic][inb]*g1; 
			complexv G1Scaaa = Mijkl.data5D[ic][ia][ia][ia][inb]*g1; 

			// g 2 derivative functions
			complexv G2Scaac = Mijkl.data5D[ic][ia][ia][ic][inb]*g2; 
            
            complexv cot = coni*(evals.data1D[ic]-evals.data1D[ia])*ti;
	
			g0 = G1Saaac-G1Sccac +2.0*coni*laccc;
			g1 = G1Scacc-G1Scaaa +2.0*coni*lcacc;
			g2 = exp( -G0Saaaa-G0Scccc + G0Sccaa+G0Saacc - 2.0*coni*ti*(lcccc-lccaa)  + cot);

			fval  +=  ( -g0*g1 + G2Scaac ) * g2 *step;
		}
        }

        rates.data2D[ia][ic] = 2.0*fval.real();
                    }



        // filling diagonal values
        for(int ia = 0; ia<num; ia++)
        {
        	rates.data2D[ia][ia] = 0.0;
        	for(int ib = 0; ib<num; ib++)
        		if(ia != ib )
        		{
        			rates.data2D[ia][ia] -= rates.data2D[ib][ia];
        		}
        }



        cout<<"# Redfield (Mod) rate matrix:\n";
        for(int ia = 0; ia<num; ia++)
        {
        	for(int ib = 0; ib<num; ib++)
        	{
        		if(ib == num-1)
        			cout<<rates.data2D[ia][ib]<<"\n";
        		else
        			cout<<rates.data2D[ia][ib]<<"\t";
        	}
        }
        
    }
    
    
    return rates;


}


storage<complexv> calculator_redfield::AddLifetimeDephasings()
{
    // here I create redfield rates and 
    // dephasing constants from mfun
    
    int numG = evals.GetSize();
        
    // number of bath oscillators
    int numBosc = mfun.GetSize();
    
    // next I make dephasings
    if( !dephasings.IsAlloc() )
        dephasings.Allocate(numG,numG);
    
    
    cout<<"# Lifetime dephasing matrix:\n";
    for(int ia = 0; ia<numG; ia++)
        for(int ib = 0; ib<numG; ib++)
        {
            complexv value = 0.0;
            //cout<<"ia ib: "<<ia<<" "<<ib<<" "<<kij.data3D[ia][ib][0]<<"\n";
            if(ia != ib )
            {
                for(int ic=0; ic<numG; ic++)
                {
                    
                    for(int iosc = 0; iosc<numBosc; iosc++)
                    {
                        // adding off diagonal fluctuations
                        if(ic != ia)
                        {
                        	double omegaij = evals.data1D[ia]-evals.data1D[ic];
                        	//omegaij += -reorganizations.data2D[ia][ia] + reorganizations.data2D[ic][ic];

                        	value
                            += kij.data3D[ia][ic][iosc]*mfun.data1D[iosc].Get(omegaij);
                            //cout<<"ia "<< ia<<" ib "<< ib<<" ic "<< ic<<"\n";
                            //cout<<"value "<< value<<"\n";
                        }
                        if(ic != ib)
                        {
                        	double omegaij = evals.data1D[ib]-evals.data1D[ic];
                        	//omegaij += -reorganizations.data2D[ib][ib] + reorganizations.data2D[ic][ic];

                        	value
                            += kij.data3D[ib][ic][iosc]*(mfun.data1D[iosc].Get(omegaij).conjugate());
                            //cout<<"ia "<< ia<<" ib "<< ib<<" ic "<< ic<<"\n";
                            //cout<<"value "<< value<<"\n";
                        }
                    }
                }
                
                dephasings.data2D[ia][ib] += value;
                
                
            }
        
            
            
            cout<<value.real()<<" "<<value.imag();
            if(ib==numG-1) cout<<"\n";
            else cout<<"\t";
            
            //cout<<"ia ib: "<<ia<<" "<<ib<<" "<<evals.data1D[ia]<<" "<<evals.data1D[ib]<<" "<<dephasings.data2D[ia][ib]<<"\n";
        }
   

    return dephasings;
}

storage<complexv> calculator_redfield::AddPureDephasings()
{
    // here I create redfield rates and 
    // dephasing constants from mfun
    
    int numG = evals.GetSize();
    
    // number of bath oscillators
    int numBosc = mfun.GetSize();
    
    
    // next I make dephasings
    if( !dephasings.IsAlloc() )
        dephasings.Allocate(numG,numG);
    
    cout<<"# Pure dephasing matrix:\n";
    for(int ia = 0; ia<numG; ia++)
        for(int ib = 0; ib<numG; ib++)
        {
            complexv value = 0.0;
            //cout<<"ia ib: "<<ia<<" "<<ib<<" "<<kij.data3D[0][ia][ib]<<"\n";
            if(ia != ib )
            {
                //for(int ic=0; ic<numG; ic++)
                //{
                    
                    for(int iosc = 0; iosc<numBosc; iosc++)
                    {
                        // adding off diagonal fluctuations
                    	complexv mfunv = mfun.data1D[iosc].Get(0);
                       
                        value 
                            += gij.data3D[ia][ia][iosc]*mfunv;
                        
                        value
                            += gij.data3D[ib][ib][iosc]*mfunv.conjugate();

                        value
                            -= 2.0*gij.data3D[ib][ia][iosc]*mfunv.real();
                    }
                //}
                
                dephasings.data2D[ia][ib] += value;
            }
            
            
            
            
            cout<<value.real()<<" "<<value.imag();
            if(ib==numG-1)
                cout<<"\n";
            else
                cout<<"\t";

        }
    
    
    return dephasings;

    
}

storage<complexv> calculator_redfield::GetRelaxationSuperoperatorM(int block)
{
    // Markovian superoperator
	// R_{ab,cd}

    if(supermatrix.IsAlloc())
    {
        supermatrix.Delete();
    }

    if(!evals.IsAlloc())
    {
        cout<<"Error: calculator_redfield::GetRelaxationSuperoperatorM(): eigenvalues are missing\n";
        return supermatrix;
    }
    
    if(!Mijkl.IsAlloc())
    {
        cout<<"Error: calculator_redfield::GetRelaxationSuperoperatorM(): fluctuation matrix is missing\n";
        return supermatrix;
    }
    
    // number of bath oscillators
    int numBosc = mfun.GetSize();

    
    
    
    
    // calculating reorganization energies from M (if necessary)
    // calculaTING REOGANIZATION functions
    if(!reorganizations.IsAlloc())
    { 
        int num = evals.GetSize();
        reorganizations.Allocate(num,num);
        cout<<"# reorganization energies:\n";

        // calculating the partition function
        for(int ind = 0; ind<num; ind++)
            for(int inr = 0; inr<num; inr++)
            {

                double part = 0;

                for(int inb = 0; inb<numBosc; inb++)
                {
                    part += Mijkl.data5D[ind][ind][inr][inr][inb]*(mfun.data1D[inb].Get(0)).imag();
                }
                reorganizations.data2D[ind][inr] = -part;

                if(inr == num-1)
                    cout<<reorganizations.data2D[ind][inr]<<"\n";
                else
                    cout<<reorganizations.data2D[ind][inr]<<"\t";
            }
    }
    
    int numL,numR,shL,shR;
    if(block == 0)
    {    
        numL = numG;
        numR = numG;
        shL = 0;
        shR = 0;
    }
    else if(block == 10)
    {    
        numL = numE;
        numR = numG;
        shL = numG;
        shR = 0;
    }
    else if(block == 11)
    {    
        numL = numE;
        numR = numE;
        shL = numG;
        shR = numG;
    }
    else if(block == 20)
    {    
        numL = numF;
        numR = numG;
        shL = numG+numE;
        shR = 0;
    }
    else if(block == 21)
    {    
        numL = numF;
        numR = numE;
        shL = numG+numE;
        shR = numG;
    }
    else if(block == 22)
    {    
        numL = numF;
        numR = numF;
        shL = numG+numE;
        shR = numG+numE;
    }
    else if(block == 123123)
    {    
        numL = numG+numE+numF;
        numR = numG+numE+numF;
        shL = 0;
        shR = 0;
    }

    
    supermatrix.Allocate(numL,numR,numL,numR);
    if (flagLindblad == 1) return GetRelaxationSuperoperatorLindblad1(block);

    for(int ia = 0; ia<numL; ia++)
        for(int ib = 0; ib<numR; ib++)
            for(int ia1 = 0; ia1<numL; ia1++)
                for(int ib1 = 0; ib1<numR; ib1++)
        {
            complexv value = 0.0;

            for(int iosc = 0; iosc<numBosc; iosc++)
            {


                if(ib == ib1)
                {
            	    for(int ic=0; ic<numL; ic++)
                    {
                       	double omegaij = evals.data1D[ia1+shL]-evals.data1D[ic+shL];
                      	//omegaij += -reorganizations.data2D[ia1+shL][ia1+shL] + reorganizations.data2D[ic+shL][ic+shL];
                       	value += Mijkl.data5D[ia+shL][ic+shL][ic+shL][ia1+shL][iosc]*mfun.data1D[iosc].Get(omegaij);
                    }
                }
                if(ia == ia1)
                {
            	    for(int ic=0; ic<numR; ic++)
                    {
                       	double omegaij = evals.data1D[ib1+shR]-evals.data1D[ic+shR];
                      	//omegaij += -reorganizations.data2D[ib1+shR][ib1+shR] + reorganizations.data2D[ic+shR][ic+shR];
                       	value += complexv(Mijkl.data5D[ib1+shR][ic+shR][ic+shR][ib+shR][iosc]*mfun.data1D[iosc].Get(omegaij)).conjugate();
                    }
                }
                if(true)
                {
                  	double omegaij = evals.data1D[ib1+shR]-evals.data1D[ib+shR];
                  	//omegaij += -reorganizations.data2D[ib1+shR][ib1+shR] + reorganizations.data2D[ib+shR][ib+shR];
                  	value -= complexv(Mijkl.data5D[ib1+shR][ib+shR][ia+shL][ia1+shL][iosc]*mfun.data1D[iosc].Get(omegaij)).conjugate();
                }
                if(true)
                {
                  	double omegaij = evals.data1D[ia1+shL]-evals.data1D[ia+shL];
                  	//omegaij += -reorganizations.data2D[ia1+shL][ia1+shL] + reorganizations.data2D[ia+shL][ia+shL];
                  	value -= Mijkl.data5D[ib1+shR][ib+shR][ia+shL][ia1+shL][iosc]*mfun.data1D[iosc].Get(omegaij);
                }
                supermatrix.data4D[ia][ib][ia1][ib1] -= value;
            }
        }


    return supermatrix;
}

storage<complexv> 
calculator_redfield::GetRelaxationSuperoperatorLindblad1(int block)
// scheme based on proper relaxation rates
{

    int numL,numR,shL,shR;
    if(block == 0)
    {    
        numL = numG;
        numR = numG;
        shL = 0;
        shR = 0;
    }
    else if(block == 10)
    {    
        numL = numE;
        numR = numG;
        shL = numG;
        shR = 0;
    }
    else if(block == 11)
    {    
        numL = numE;
        numR = numE;
        shL = numG;
        shR = numG;
    }
    else if(block == 20)
    {    
        numL = numF;
        numR = numG;
        shL = numG+numE;
        shR = 0;
    }
    else if(block == 21)
    {    
        numL = numF;
        numR = numE;
        shL = numG+numE;
        shR = numG;
    }
    else if(block == 22)
    {    
        numL = numF;
        numR = numF;
        shL = numG+numE;
        shR = numG+numE;
    }
    else if(block == 123123)
    {    
        numL = numG+numE+numF;
        numR = numG+numE+numF;
        shL = 0;
        shR = 0;
    }



    // number of bath oscillators
    int numBosc = mfun.GetSize();

    int num = numG+numE+numF;
    
	// calculating supermatrix Zabcd - the Lindblad correlation supermatrix
        // from Redfield-like relations
    storage<complexv> zabcd(4);
    zabcd.Allocate(num,num,num,num);

    storage<complexv> zV(3);
    zV.Allocate(num,num,numBosc);

    for(int ia = 0; ia<num; ia++)
    for(int ic = 0; ic<num; ic++)
    for(int iosc = 0; iosc<numBosc; iosc++){
		zV.data3D[ia][ic][iosc] = sqrt(Mijkl.data5D[ic][ia][ia][ic][iosc]);
    }

// calulcating autocorrelations
    for(int ib = 0; ib<num; ib++)
    for(int id = 0; id<num; id++)
    for(int ia = 0; ia<num; ia++)
    for(int ic = 0; ic<num; ic++)
    for(int iosc = 0; iosc<numBosc; iosc++){

        double omegaij = evals.data1D[id]-evals.data1D[ib];
// 		//omegaij = evals.data1D[id]-evals.data1D[ib] -reorganizations.data2D[id][id] + reorganizations.data2D[ib][ib];
        complexv value1 = mfun.data1D[iosc].Get(omegaij);
                
 		omegaij = evals.data1D[ic]-evals.data1D[ia];
// 		//omegaij = evals.data1D[ic]-evals.data1D[ia] -reorganizations.data2D[ic][ic] + reorganizations.data2D[ia][ia];
 		complexv value2 = mfun.data1D[iosc].Get(omegaij);


            zabcd.data4D[ib][id][ia][ic] +=  (zV.data3D[ib][id][iosc]).conjugate()* zV.data3D[ia][ic][iosc] * (value1.conjugate()+value2);
			cout<<"autocorrelation: "<<ib<<" "<<id<<" "<<ia<<" "<<ic<<" "<<zabcd.data4D[ib][id][ia][ic]<<"\n";
        }

// // calulcating autocorrelations
//     for(int ib = 0; ib<num; ib++)
//     for(int id = 0; id<num; id++)
//     for(int ia = 0; ia<num; ia++)
//     for(int ic = 0; ic<num; ic++)
//         {
//        		int ib = ia;
//  		int id = ic;
//             complexv value1 = 0.0;
//             complexv value2 = 0.0;
// 		double omegaij;

//             for(int iosc = 0; iosc<numBosc; iosc++)
//             {
// 		omegaij = evals.data1D[id]-evals.data1D[ib];
// 		//omegaij = evals.data1D[id]-evals.data1D[ib] -reorganizations.data2D[id][id] + reorganizations.data2D[ib][ib];
// 		value1 += Mijkl.data5D[id][ib][ia][ic][iosc]*mfun.data1D[iosc].Get(omegaij);
                
// 		omegaij = evals.data1D[ic]-evals.data1D[ia];
// 		//omegaij = evals.data1D[ic]-evals.data1D[ia] -reorganizations.data2D[ic][ic] + reorganizations.data2D[ia][ia];
// 		value2 += Mijkl.data5D[id][ib][ia][ic][iosc]*mfun.data1D[iosc].Get(omegaij);
// 	    }
//             zabcd.data4D[ib][id][ia][ic] = value1.conjugate()+value2;


// 			cout<<"autocorrelation: "<<ib<<" "<<id<<" "<<ia<<" "<<ic<<" "<<zabcd.data4D[ib][id][ia][ic]<<"\n";

//         }


//  // calulcating correlations
//    for(int ia = 0; ia<num; ia++)
//        for(int ib = 0; ib<num; ib++)
//            for(int ic = 0; ic<num; ic++)
//                for(int id = 0; id<num; id++)
//      		if(ib != ia ||  id != ic)
//		{
//	            zabcd.data4D[ib][id][ia][ic] = sqrt(zabcd.data4D[ib][id][ib][id]*zabcd.data4D[ia][ic][ia][ic]);
//		}

    
	// calculating relaxation supermatrix from the Lindblad correlation supermatrix
    for(int ia = 0; ia<numL; ia++)
        for(int ib = 0; ib<numR; ib++)
            for(int ic = 0; ic<numL; ic++)
                for(int id = 0; id<numR; id++)
        {
            complexv value = 0.0;

                    if(ib == id)                    
                	for(int ie=0; ie<numL; ie++)
			value = value - 0.5*zabcd.data4D[ie+shL][ia+shL][ie+shL][ic+shL];

                    if(ia == ic)
                	for(int ie=0; ie<numR; ie++)
			value = value - 0.5*zabcd.data4D[ie+shR][id+shR][ie+shR][ib+shR];

                supermatrix.data4D[ia][ib][ic][id] = (value + zabcd.data4D[ib+shR][id+shR][ia+shL][ic+shL]);
        }


    return supermatrix;
}

//void calculator_redfield::AddMijkl(storage<complexv>& iM)
//{
//	if(iM.CheckDimension() == 5)
//		Mijkl = iM;
//	else
//	{
//		cout<<"Error: calculator_redfield::AddMijkl wrong dimensionality of M\n";
//	}
//
//}
void calculator_redfield::AddMijkl(storage<double>& iM)
{
	if(iM.CheckDimension() == 5)
    {
		Mijkl = iM;
	}
    else
	{
		cout<<"Error: calculator_redfield::AddMijkl wrong dimensionality of M\n";
	}

}

// Schrodinger picture Redfield memory kernel
storage<complexv> calculator_redfield::GetMemoryKernel(double time, double& deltat, int& numT, int block)
// block should be 
// 0, 10, 11, 20, 21, 22
{

    
//    if(!kernel.IsAlloc())
//    {
        //cout<<"Error: calculator_redfield::GetMemoryKernelInteraction: kernel must be allocated\n";
        //return kernel;
//    }

    if(!evals.IsAlloc())
    {
        cout<<"Error: calculator_redfield::SetupMemoryKernel: eigenvalues are missing\n";
        return kernel;
    }
    
    if(!Mijkl.IsAlloc())
    {
        cout<<"Error: calculator_redfield::SetupMemoryKernel: fluctuation matrix is missing\n";
        return kernel;
    }


    if(numT == 0)
	{
		// define history parameters from the zeroth correlation function
		deltat = cfun.data1D[0].GetStep();
		numT = cfun.data1D[0].GetN();
	}

    int numosc = cfun.GetSize();




    int numL,numR,shL,shR;
    if(block == 0)
    {    
        numL = numG;
        numR = numG;
        shL = 0;
        shR = 0;
    }
    else if(block == 10)
    {    
        numL = numE;
        numR = numG;
        shL = numG;
        shR = 0;
    }
    else if(block == 11)
    {    
        numL = numE;
        numR = numE;
        shL = numG;
        shR = numG;
    }
    else if(block == 20)
    {    
        numL = numF;
        numR = numG;
        shL = numG+numE;
        shR = 0;
    }
    else if(block == 21)
    {    
        numL = numF;
        numR = numE;
        shL = numG+numE;
        shR = numG;
    }
    else if(block == 22)
    {    
        numL = numF;
        numR = numF;
        shL = numG+numE;
        shR = numG+numE;
    }
    else if(block == 123123)
    {    
        numL = numG+numE+numF;
        numR = numG+numE+numF;
        shL = 0;
        shR = 0;
    }


    cout<<"numL = "<<numL;
    cout<<"shL = "<<shL;
    cout<<"numR = "<<numR;
    cout<<"shR = "<<shR;


      if(!kernel.IsAlloc())
          kernel.Allocate(numT,numL, numR, numL, numR);


      for(int it=0;it<numT;it++)
      {
        double tau = it*deltat;

        for(int ia=0;ia<numL;ia++)
        for(int ia1=0;ia1<numL;ia1++)

        for(int ib=0;ib<numR;ib++)
        for(int ib1=0;ib1<numR;ib1++)
	{
		complexv R1=0;
		if(ib==ib1)
	        for(int ic=0;ic<numL;ic++)
		{
			double omegat = evals.data1D[ic+shL]-evals.data1D[ib1+shR];
		        for(int io=0;io<numosc;io++)
			R1 += exp(cnni*omegat*tau)*Mijkl.data5D[ia+shL][ic+shL][ic+shL][ia1+shL][io]*cfun.data1D[io].Get(tau);
		}
		complexv R4=0;
		if(ia==ia1)
	        for(int id=0;id<numR;id++)
		{
			double omegat = evals.data1D[ia1+shL]-evals.data1D[id+shR];
		        for(int io=0;io<numosc;io++)
			R4 += exp(cnni*omegat*tau)*Mijkl.data5D[ib1+shR][id+shR][id+shR][ib+shR][io]*conj(cfun.data1D[io].Get(tau));
		}
		complexv R2=0;
		{
			double omegat = evals.data1D[ia1+shL]-evals.data1D[ib+shR];
		        for(int io=0;io<numosc;io++)
			R2 += exp(cnni*omegat*tau)*Mijkl.data5D[ib1+shR][ib+shR][ia+shL][ia1+shL][io]*conj(cfun.data1D[io].Get(tau));
		}
		complexv R3=0;
		{
			double omegat = evals.data1D[ia+shL]-evals.data1D[ib1+shR];
		        for(int io=0;io<numosc;io++)
			R3 += exp(cnni*omegat*tau)*Mijkl.data5D[ia+shL][ia1+shL][ib1+shR][ib+shR][io]*cfun.data1D[io].Get(tau);
		}
	
        if(it==0)
            kernel.data5D[it][ia][ib][ia1][ib1] = -(R1-R2-R3+R4)/2.0;
            // this is needed for integration by trapecia
        else
            kernel.data5D[it][ia][ib][ia1][ib1] = -(R1-R2-R3+R4);
      }
    }

//    cout<<"M(omega=0) = "<<mfun.data1D[0].Get(0)<<"\n";
//    cout<<"M(t) function:\n";
//    cout<<"dt = "<<deltat<<"\n";
//    for(int it=0;it<numT;it++)
//    cout<<cfun.data1D[0].Get(deltat*it)<<"\n";

    cout<<"Kernel specific element kernel.data5D[it][1][1][2][2]\n";
    cout<<"Mijkl.data5D[1][2][2][1][0] = "<<Mijkl.data5D[1+1][2+1][2+1][1+1][0]<<"\n";
    cout<<"dt = "<<deltat<<"\n";
    for(int it=0;it<numT;it++)
    cout<<kernel.data5D[it][1][1][2][2]<<"\n";
        
    return kernel;
            
}

////////////////////////////
// the factorized form of the memory kernel (without cfun)
storage<complexv> calculator_redfield::GetMemoryKern(int block)
// block should be 
// 0, 10, 11, 20, 21, 22
{

    if(!Mijkl.IsAlloc())
    {
        cout<<"Error: calculator_redfield::SetupMemoryKernel: fluctuation matrix is missing\n";
        return kernel;
    }



    int numL,numR,shL,shR;
    if(block == 0)
    {    
        numL = numG;
        numR = numG;
        shL = 0;
        shR = 0;
    }
    else if( (block == 10) || (block == 11) )
    {    
        numL = numE+numG;
        numR = numE+numG;
        shL = 0;
        shR = 0;
    }
    else 
    {    
        numL = numF+numE+numG;
        numR = numF+numE+numG;
        shL = 0;
        shR = 0;
    }

    // here M matrix is converted into blocks and returned
    storage<complexv> retM(3);
    retM.Allocate(cfun.GetSize(),numL*numL,numR*numR);

    for(int ia=0;ia<numL;ia++)
    for(int ia1=0;ia1<numL;ia1++)

    for(int ib=0;ib<numR;ib++)
    for(int ib1=0;ib1<numR;ib1++)

    for(int io=0;io<cfun.GetSize();io++)
        retM.data3D[io][ia*numL+ia1][ib*numR+ib1] =  Mijkl.data5D[ia][ia1][ib][ib1][io];
        
    return retM;
            
}


//------------------------------------

void calculator_redfield::Setup()
{
    // general flag
    ready = 0;

    // system parameters
    evals.SetDimension(1);

    rates.SetDimension(2);    
    dephasings.SetDimension(2);
    reorganizations.SetDimension(2);
    
    mfun.SetDimension(1); // for a set of spectral densities 
    cfun.SetDimension(1); // for a set of spectral densities 
    gij.SetDimension(3); // amplitudes for all pairs
    kij.SetDimension(3); // amplitudes for all pairs
    Mijkl.SetDimension(5); // Hilbert notation is used

    supermatrix.SetDimension(4); // Hilbert notation is used
    kernel.SetDimension(5); // 
    // transport rates:

    flagLindblad = 0;
	tempr = 0;
	//nonsecular = 0;
}


