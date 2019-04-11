// everything is performed in terms of the accessory bath functions a1,a2 and a3
// at the moment that algorithm is most advanced
#include"../storage/storage.hpp"
#include"eigen.hpp"
#include"include-LAPACK/dsyev.h"
#include"include-LAPACK/dgeev.h"
#include"include-LAPACK/zgeev.h"
#include"../constants/constants.hpp"

eigen::eigen( int setup)
{
    //ltype = setup;
}
eigen::eigen()
{
    //ltype = 0;
    
}
eigen::~eigen()
{
    //ltype = 0;
    
}

// real symmetric matrix
int eigen::GetEigenSolution(storage<double>& imatrix, storage<double>& evals, storage<double>& evecs)
{
    return GetEigenSolution(&imatrix, &evals, &evecs);
}

// real symmetric matrix
int eigen::GetEigenSolution(storage<double>* imatrix, storage<double>* evals, storage<double>* evecs)
{
    
    // getting the number of elements
    int Num;
    evals->GetSize(Num);
    if(Num == 0)
    {
        cout<<"Warning: eigenvalue solution for zero size matrix? ... skipping\n";
        return 1;
    }

    cout<<"# Calculating eigenvalues\n";

    double* Evecs = new double[Num*Num];

	for(int ii=0;ii<Num;ii++)
        for(int ij=0;ij<Num;ij++)
            Evecs[ii*Num+ij]=imatrix->data2D[ii][ij];
    
    int LWORK;
    double *WORK;
    int INFO;
    WORK = new double[2];
	LWORK=-1;
	dsyev_( 'V', 'U',Num,Evecs,Num,evals->data1D,WORK,LWORK,INFO );
	LWORK= (int) WORK[0];
	delete[] WORK;
	WORK = new double[LWORK];
	dsyev_( 'V', 'U',Num,Evecs,Num,evals->data1D,WORK,LWORK,INFO );
    delete[] WORK;
    
    if(INFO != 0)
    {
        cout<<"Error: eigenvalue solution failed ... skipping\n";
        
        delete[] Evecs;
        return INFO;
    }

	
    cout<<"# Eigenvalues done.\n";
// Eigenvalues:
//***    cout<<"# Values:\n";
//***    for(int ii=0;ii<Num;ii++)
//***       cout<<evals->data1D[ii]<<"\n";
    
    
    
    for(int ii=0;ii<Num;ii++)
        for(int ij=0;ij<Num;ij++)
        {
            //if(ltype)
                evecs->data2D[ii][ij]=Evecs[ij*Num+ii];
            //else
            //    evecs->data2D[ii][ij]=Evecs[ii*Num+ij];
        }
    delete[] Evecs;


    
    //*** cout<<"# Orthogonality:\n";
    // renormalization
    // experience is that sometimes the eigenvectors are not normalized to 1:
   for(int ii=0;ii<Num;ii++)
    {
            double nor = 0.0;

            // calculating norm and rescaling
                for(int ik=0;ik<Num;ik++)
                    nor += evecs->data2D[ik][ii]*evecs->data2D[ik][ii];

                for(int ik=0;ik<Num;ik++)
                {
                    evecs->data2D[ik][ii] /=nor;
                }
   }
    //***for(int ii=0;ii<Num;ii++)
    //***    for(int ij=0;ij<Num;ij++)
    //***    {
    //***        double nor = 0.0;
    //***            // checking orthogonality
    //***            for(int ik=0;ik<Num;ik++)
    //***                nor += evecs->data2D[ik][ii]*evecs->data2D[ik][ij];
    //***
    //***        cout<<nor;
    //***        if(ij == Num-1) cout<<"\n";
    //***        else cout<<"\t";
    //***    }
    
    
    
    



    
    return INFO;

}	


// square non-symmetric real matrix imatrix
int eigen::GetEigenSolution(storage<double>& imatrix, storage<complexv>& evals, storage<complexv>& evecsL,  storage<complexv>& evecsR)
{

    // getting the number of elements
    int Num;
    evals.GetSize(Num);
    if(Num == 0)
    {
        cout<<"Warning: eigenvalue solution for zero size matrix? ... skipping\n";
        return 1;
    }

    cout<<"# Calculating eigenvalues\n";

    double* matrix = new double[Num*Num];

	for(int ii=0;ii<Num;ii++)
        for(int ij=0;ij<Num;ij++)
	{
                //if(ltype)
            		matrix[ii*Num+ij]=imatrix.data2D[ij][ii];
		//else
            	//	matrix[ii*Num+ij]=imatrix.data2D[ii][ij];
	}

    double* valsR = new double[Num];
    double* valsI = new double[Num];

    double* vL = new double[Num*Num];
    double* vR = new double[Num*Num];

	int iNFO;
	int iwORK = -1;
	double dwORK;
	dgeev_( 'V','V',Num,matrix, Num, valsR, valsI, vL,Num, vR, Num,  &dwORK,iwORK,iNFO );

	iwORK= (int) dwORK;
    double* wORK = new double[iwORK];
	dgeev_( 'V','V',Num,matrix, Num, valsR, valsI, vL,Num, vR, Num,  wORK,iwORK,iNFO );
    delete[] wORK;
    delete[] matrix;
    
    if(iNFO != 0)
    {
        cout<<"Error: eigenvalue solution failed ... skipping\n";
        
        delete[] vR;
        delete[] vL;
        delete[] valsI;
        delete[] valsR;

        return iNFO;
    }


    cout<<"# Eigenvalues done.\n";

    for(int ii=0;ii<Num;ii++)
    {
    	evals.data1D[ii] = complexv(valsR[ii],valsI[ii]);

    	if(fabs(valsI[ii])<constants::smallepsilon)
    	{
    		// eigenvalue is real
            for(int ij=0;ij<Num;ij++)
            {
                //if(ltype)
                    evecsR.data2D[ii][ij]=vR[ij*Num+ii];
                //else
                //    evecsR.data2D[ii][ij]=vR[ii*Num+ij];

                //if(ltype)
                    evecsL.data2D[ii][ij]=vL[ij*Num+ii];
                //else
                //    evecsL.data2D[ii][ij]=vL[ii*Num+ij];
            }
    	}
    	else
    	{
    		// eigenvalue is complex
            evals.data1D[ii+1] = complexv(valsR[ii+1],valsI[ii+1]);

    		// adding two eigenvectors
            for(int ij=0;ij<Num;ij++)
            {
                // this is written only for ltype == 1
            		evecsR.data2D[ij][ii]=complexv(vR[ii*Num+ij],vR[(ii+1)*Num+ij]);
            		evecsR.data2D[ij][ii+1]=complexv(vR[ii*Num+ij],-vR[(ii+1)*Num+ij]);


            		evecsL.data2D[ij][ii]=complexv(vL[ii*Num+ij],vL[(ii+1)*Num+ij]);
            		evecsL.data2D[ij][ii+1]=complexv(vL[ii*Num+ij],-vL[(ii+1)*Num+ij]);
            }
            ii++;
    	}
    }
    
    // Eigenvalues:
//***    cout<<"# Values:\n";
//***    for(int ii=0;ii<Num;ii++)
//***       cout<<evals.data1D[ii]<<"\n";


//    cout<<"# Eigenvalue equation:\n";
//    for(int ii=0;ii<Num;ii++)
//        for(int ij=0;ij<Num;ij++)
//        {
//
//		cout<<ii<<" "<<ij<<" Au: ";
//
//            complexv nor = 0.0;
//            {
//                // checking orthogonality
//                for(int ik=0;ik<Num;ik++)
//                    nor += imatrix.data2D[ii][ik]*evecsR.data2D[ik][ij];
//            }
//            cout<<nor;
//
//           nor =  evecsR.data2D[ii][ij]*evals.data1D[ij];
//           
//	cout<<" ul: ";
//            cout<<nor;
//
//	cout<<"\n";
//
//        }
    

    
   //*** cout<<"# Orthogonality:\n";


    
    // renormalization
    // experience is that sometimes the eigenvectors are not normalized to 1:
        for(int ii=0;ii<Num;ii++)
        {
            complexv nor = 0.0;
                // calculating norm and rescaling
                for(int ik=0;ik<Num;ik++)
                    nor += evecsL.data2D[ik][ii]*evecsR.data2D[ik][ii];
                nor = nor.abs();
                nor = sqrt(nor);
                for(int ik=0;ik<Num;ik++)
                {
                    evecsL.data2D[ik][ii] /=nor;
                    evecsR.data2D[ik][ii] /=nor;
                }
        }
//***    for(int ii=0;ii<Num;ii++)
//***        for(int ij=0;ij<Num;ij++)
//***        {
//***            complexv nor = 0.0;
//***                // checking orthogonality
//***                for(int ik=0;ik<Num;ik++)
//***                    nor += evecsL.data2D[ik][ii]*evecsR.data2D[ik][ij];
//***            cout<<nor;
//***            if(ij == Num-1) cout<<"\n";
//***            else cout<<"\t";
//***        }
    
    
    
    
    
    
    delete[] vR;
    delete[] vL;
    delete[] valsI;
    delete[] valsR;

    return iNFO;

}


// square non-symmetric complex matrix imatrix
int eigen::GetEigenSolution(storage<complexv>& imatrix, storage<complexv>& evals, storage<complexv>& evecsL,  storage<complexv>& evecsR)
{

    // getting the number of elements
    int Num;
    evals.GetSize(Num);
    if(Num == 0)
    {
        cout<<"Warning: eigenvalue solution for zero size matrix? ... skipping\n";
        return 1;
    }

    cout<<"# Calculating eigenvalues\n";

    std::complex<double>* matrixA = new std::complex<double>[Num*Num];
    //std::complex<double>* matrixB = new std::complex<double>[Num*Num];

	for(int ii=0;ii<Num;ii++)
        for(int ij=0;ij<Num;ij++)
	{
                //if(ltype)
                //{
            		matrixA[ii*Num+ij]=imatrix.data2D[ij][ii];
            		//matrixB[ii*Num+ij]=imatrix.data2D[ij][ii];
                //}
		//else
                //{
            	//	matrixA[ii*Num+ij]=imatrix.data2D[ii][ij];
            	//	//matrixB[ii*Num+ij]=imatrix.data2D[ii][ij];
                //}
	}

    std::complex<double>* valsA = new std::complex<double>[Num];
    //std::complex<double>* valsB = new std::complex<double>[Num];

    std::complex<double>* vL = new std::complex<double>[Num*Num];
    std::complex<double>* vR = new std::complex<double>[Num*Num];

	int iNFO;
	int iwORK = -1;
	std::complex<double> qwORK;
	double* dwORK = new double[2*Num];
	zgeev_( 'V','V',Num,matrixA,Num,valsA, vL,Num, vR, Num, &qwORK,iwORK,dwORK, iNFO );

	iwORK= (int) qwORK.real();
        std::complex<double>* cwORK = new std::complex<double>[iwORK];
	zgeev_( 'V','V',Num,matrixA,Num,valsA, vL,Num, vR, Num, cwORK,iwORK,dwORK, iNFO );
    delete[] dwORK;
    delete[] cwORK;
    delete[] matrixA;
    //delete[] matrixB;

    
    
    if(iNFO != 0)
    {
        cout<<"Error: eigenvalue solution failed ... skipping\n";
        
        delete[] vR;
        delete[] vL;
        delete[] valsA;
        //delete[] valsB;
        
        return iNFO;
    }
    
    
    cout<<"# Eigenvalues done.\n";

    for(int ii=0;ii<Num;ii++)
    {
        //if( std::abs(valsA[ii])<constants::smallepsilon)
        //    evals.data1D[ii] = 0.0;
        //else if(std::abs(valsB[ii])<constants::smallepsilon)
        //{
        //    cout<<"Error: eigenvalue solution generated diverging eigenvalues ... result is not controlled\n";
        //    evals.data1D[ii] = 1.0/constants::smallepsilon;
        //}
        //        else
            evals.data1D[ii] = valsA[ii];


        
            for(int ij=0;ij<Num;ij++)
            {
                //if(ltype)
                    evecsR.data2D[ii][ij]=vR[ij*Num+ii];
                //else
                //    evecsR.data2D[ii][ij]=vR[ii*Num+ij];

                //if(ltype)
                    evecsL.data2D[ii][ij]=conj(vL[ij*Num+ii]);
                //else
                //    evecsL.data2D[ii][ij]=vL[ii*Num+ij];
            }
    }
    
        // Eigenvalues:
//***    cout<<"# Values:\n";
//***    for(int ii=0;ii<Num;ii++)
//***       cout<<evals.data1D[ii]<<"\n";
    
    //***cout<<"# Orthogonality:\n";

    
    
    // renormalization
    // experience is that sometimes the eigenvectors are not normalized to 1:
    for(int ii=0;ii<Num;ii++)
        {
            complexv nor = 0.0;
                // calculating norm and rescaling
                for(int ik=0;ik<Num;ik++)
                    nor += evecsL.data2D[ik][ii]*evecsR.data2D[ik][ii];
                //nor = nor.abs();
                nor = sqrt(nor);
                for(int ik=0;ik<Num;ik++)
                {
                    evecsL.data2D[ik][ii] /=nor;
                    evecsR.data2D[ik][ii] /=nor;
                }
    }
//    for(int ii=0;ii<Num;ii++)
//        for(int ij=0;ij<Num;ij++)
//        {
//            complexv nor = 0.0;
//                // checking orthogonality
//                for(int ik=0;ik<Num;ik++)
//                    nor += evecsL.data2D[ik][ii]*evecsR.data2D[ik][ij];
//                
//            cout<<nor;
//            if(ij == Num-1) cout<<"\n";
//            else cout<<"\t";
//        }
    
    
    
    
    
    delete[] vR;
    delete[] vL;
        delete[] valsA;
        //delete[] valsB;

    return iNFO;

}
