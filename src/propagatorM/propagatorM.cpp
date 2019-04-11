#include "propagatorM.hpp"
#include"../toolsInterpolate/toolsInterpolate.hpp"
#include "../storage/storage.hpp"
#include "../constants/constants.hpp"
#include "../eigen/eigen.hpp"
#include<cmath>


// static members
 //int propagatorM::dimension;
 //int* propagatorM::dataInt;
 //double* propagatorM::dataset1D;
 //double** propagatorM::dataset2D;


 propagatorM::propagatorM()
 {
	 dimension = 0;
	 dataInt = 0;
	 //dataset1D.SetDimension(1);
	 mastermatrix.SetDimension(2);
	 mastermatrixComplex.SetDimension(2);

	 evals.SetDimension(1);
	 evecsR.SetDimension(2);
	 evecsL.SetDimension(2);

	configurationinitial.SetDimension(1);
	configurationinitialcomplex.SetDimension(1);
	statefinal = -1;

     pevals = 0;
     pevecsR = 0;
     pevecsL = 0;

	errcode = 0;
 }
 propagatorM::~propagatorM()
 {
	 if (dataInt != 0)
		 delete[] dataInt;
 }

// ***********************
// simple master equation
// usage 
// Init, Propagate and Exit
// int propagatorM::RightSideMaster(
//    double t, double* y, double* dy, int& num)
//{
//    for( int indl = 0; indl < num; indl ++)
//    {
//        dy[indl] = 0.0;
//        for( int indr = 0; indr < num; indr ++)
//            dy[indl] += mastermatrix.data2D[indl][indr]*y[indr];
//    }
//    return 0;
//}
//int propagatorM::InitMaster(int num, double** matrix)
//{
//    dimension = num;
//    mastermatrix.Allocate(num,num);
//    for(int i1=0; i1<num;i1++)
//        for(int i2=0; i2<num;i2++)
//    mastermatrix.data2D[i1][i2] = matrix[i1][i2];
//    dataset1D.Allocate(num*5);
//    return 0;
//}
int propagatorM::InitMaster(storage<double>& matrix)
{
    matrix.GetSize(dimension,dimension);
    mastermatrix = matrix;
    evals.Delete();
    evecsR.Delete();
    evecsL.Delete();

     pevals = 0;
     pevecsR = 0;
     pevecsL = 0;
    return 0;
}
int propagatorM::InitMaster(storage<complexv>& matrix)
{
    matrix.GetSize(dimension,dimension);
    mastermatrixComplex = matrix;
    evals.Delete();
    evecsR.Delete();
    evecsL.Delete();

     pevals = 0;
     pevecsR = 0;
     pevecsL = 0;
    return 0;
}
//int propagatorM::PropagateMaster(double* y,int numt,double dt)
//{
//    for(int ind = 0; ind< numt; ind++)
//    {
//        odeStep(ind*dt,    // current time befor the next step
//               y,   // initial values,  on output replaced by new values
//               dt,    // requested step [time dimension]
//               //(int (*)(double, double*, double*, int&))&propagatorM::RightSideMaster, // derivative function
//               dataset1D.data1D,   // for temporary data storage: size N*5
//               dimension          // number of variables
//            );
//    }
//
//    return 0;
//
//}
//int propagatorM::PropagateMaster(storage<double>& yi,int numt,double dt)
//{
//        return PropagateMaster(yi.data1D,numt,dt);
//}
int propagatorM::ExitMaster()
{
    //dataset1D.Delete();
    mastermatrix.Delete();
    mastermatrixComplex.Delete();
    evals.Delete();
    evecsR.Delete();
    evecsL.Delete();
    dimension = 0;

    pevals = 0;
    pevecsR = 0;
    pevecsL = 0;

    return 0;
}

void propagatorM::prepareGetAt(storage<double>& initials,int retindex)
{
	if(initials.CheckDimension()!=1)
	{
		cout<<"Error: propagatorM: wrong dimenision for initial conditions\n";
		errcode = 2;
	}
	else
	{
		configurationinitial = initials;
		statefinal = retindex;
	}
}
void propagatorM::prepareGetAt(storage<complexv>& initials,int retindex)
{
	if(initials.CheckDimension()!=1)
	{
		cout<<"Error: propagatorM: wrong dimenision for initial conditions\n";
		errcode = 2;
	}
	else
	{
		configurationinitialcomplex = initials;
		statefinal = retindex;
	}
}
double propagatorM::GetAt(double& time)
{
	double retvals = 0.0;

	if(errcode)
	{
		return retvals;
	}

	if(statefinal == -1)
	{
		cout<<"Error: propagatorM: set the final state\n";
		return retvals;
	}

	if ( pevecsR == 0)
	{
		 // perform eigenvalue solution
		 evals.Allocate(dimension);
		 evecsR.Allocate(dimension,dimension);
		 evecsL.Allocate(dimension,dimension);

                 // shifting the diagonal of the matrix for better orthogonality
            for(int ind = 0; ind< dimension; ind ++)
            {
                mastermatrix.data2D[ind][ind] -=1.0;
            }

	    eigen esolver(1);
	    esolver.GetEigenSolution(mastermatrix, evals, evecsL,  evecsR);
                 
            // shifting the eigenvalues back
            cout<<"# Decay constants:\n";
            for(int ind = 0; ind< dimension; ind ++)
            {
                evals.data1D[ind] +=1.0;
                cout<<evals.data1D[ind]<<"\n";
            }
                 
        
//        cout<<"evecsR:\n";
//        for(int ind = 0; ind< dimension; ind ++)
//        for(int ins = 0; ins< dimension; ins ++)
//        {
//            cout<<evecsR.data2D[ind][ins];
//            if(ins == dimension-1)
//                cout<<"\n";
//            else
//                cout<<"\t";
//        }
//        cout<<"evecsL:\n";
//        for(int ind = 0; ind< dimension; ind ++)
//            for(int ins = 0; ins< dimension; ins ++)
//            {
//                cout<<evecsL.data2D[ind][ins];
//                if(ins == dimension-1)
//                    cout<<"\n";
//                else
//                    cout<<"\t";
//            }

		// finding largest eigenvalue
		// for proper setup it must be equal to zero
		double max = -1e64;
		int maxindex;
		for(int ind = 0; ind< dimension; ind ++)
		{
			if(max < evals.data1D[ind].real())
			{
				max = evals.data1D[ind].real();
				maxindex = ind;
			}
		}
		if(fabs(max) > constants::smallepsilon)
		{
			cout<<"Error: propagatorM eigenvalue solution failed.\n";
			errcode = 1;
			return retvals;
		}
		else
		{
			// setting strong zero
			evals.data1D[maxindex] = 0.0;
		}

             pevals = evals.data1D;
             pevecsR = evecsR.data2D;
             pevecsL = evecsL.data2D;
	 }

	 // perform exponentiation and propagation
	//for(int i3=0;i3<dimension; i3++)
	//{
	int& i3 = statefinal;
		retvals = 0.0;
		 for(int i2=0;i2<dimension; i2++)
			 for(int i1=0;i1<dimension; i1++)
				 retvals += (pevecsR[i3][i2]*exp(pevals[i2]*time)*pevecsL[i1][i2]).real()*configurationinitial.data1D[i1];
	//}

	 return retvals;
}




// on input y are initial values,
// here it is accepted that propagation is real-valued
storage<double> propagatorM:: Get(storage<double>& y,double& time)
{
    if ( pevecsR != 0 )
        dimension = evals.GetSize();
    
    storage<double> retvals(1);
    retvals.Allocate(dimension);
    
	if(errcode)
		return retvals;

    // here work with pointers helps in speed
     if ( pevecsR == 0 )
     {
             // perform eigenvalue solution
             evals.Allocate(dimension);
             evecsR.Allocate(dimension,dimension);
             evecsL.Allocate(dimension,dimension);

             //int size = y.GetSize();
        for(int ind = 0; ind< dimension; ind ++)
        {
            mastermatrix.data2D[ind][ind] -=1.0;
        }

		 eigen esolver(1);
		 esolver.GetEigenSolution(mastermatrix, evals, evecsL,  evecsR);
                 
                 // shifting the eigenvalues back
                 cout<<"# Decay constants:\n";
        for(int ind = 0; ind< dimension; ind ++)
        {
            evals.data1D[ind] +=1.0;
            cout<<evals.data1D[ind]<<"\n";
        }

             // finding largest eigenvalue
             // for proper setup it must be equal to zero
             double max = -1e64;
             int maxindex;
             for(int ind = 0; ind< dimension; ind ++)
             {
			if(max < evals.data1D[ind].real())
			{
				max = evals.data1D[ind].real();
				maxindex = ind;
			}
             }
             if(fabs(max) > constants::smallepsilon)
             {
                 cout<<"Error: propagatorM eigenvalue solution failed.\n";
                 errcode = 1;
                 return retvals;
             }
             else
             {
                 // setting strong zero
                 evals.data1D[maxindex] = 0.0;
             }

             pevals = evals.data1D;
             pevecsR = evecsR.data2D;
             pevecsL = evecsL.data2D;
     }

	 // perform exponentiation and propagation
	 for(int i3=0;i3<dimension; i3++)
	{
		retvals.data1D[i3] = 0.0;
		 for(int i2=0;i2<dimension; i2++)
			 for(int i1=0;i1<dimension; i1++)
				 retvals.data1D[i3] += (pevecsR[i3][i2]*exp(pevals[i2]*time)*pevecsL[i1][i2]).real()*y.data1D[i1];
	}

	 return retvals;
}


///////////////////
// complex values
// on input y are initial values,
// here it is accepted that propagation is complex-valued
storage<complexv> propagatorM:: Get(storage<complexv>& y,double& time)
{
    if ( pevecsR != 0 )
        dimension = evals.GetSize();

    storage<complexv> retvals(1);
    retvals.Allocate(dimension);

	if(errcode)
		return retvals;

    if ( pevecsR == 0 )
    {
          // perform eigenvalue solution
          evals.Allocate(dimension);
          evecsR.Allocate(dimension,dimension);
          evecsL.Allocate(dimension,dimension);

            
          for(int ind = 0; ind< dimension; ind ++)
          {
              mastermatrixComplex.data2D[ind][ind] -=1.0;
          }

	  eigen esolver(1);
	  esolver.GetEigenSolution(mastermatrixComplex, evals, evecsL,  evecsR);
                 
          // shifting the eigenvalues back
          cout<<"# Decay constants:\n";
          for(int ind = 0; ind< dimension; ind ++)
          {
              evals.data1D[ind] +=1.0;
              cout<<evals.data1D[ind]<<"\n";
          }
//            eigen esolver(1);
//            esolver.GetEigenSolution(mastermatrixComplex, evals, evecsL,  evecsR);

        pevals = evals.data1D;
        pevecsL = evecsL.data2D;
        pevecsR = evecsR.data2D;
    }

	// perform exponentiation and propagation
	for(int i3=0;i3<dimension; i3++)
	{
		retvals.data1D[i3] = 0.0;
		 for(int i2=0;i2<dimension; i2++)
			 for(int i1=0;i1<dimension; i1++)
				 retvals.data1D[i3] += pevecsR[i3][i2]*exp(pevals[i2]*time)*pevecsL[i1][i2]*y.data1D[i1];
	}

	 return retvals;
}

complexv propagatorM::GetAtC(double& time)
{
	complexv retvals = 0.0;

	if(errcode)
	{
		return retvals;
	}

	if(statefinal == -1)
	{
		cout<<"Error: propagatorM: set the final state\n";
		return retvals;
	}

	if ( pevecsR==0 )
	{
		 // perform eigenvalue solution
		 evals.Allocate(dimension);
		 evecsR.Allocate(dimension,dimension);
		 evecsL.Allocate(dimension,dimension);

                 // shifting the diagonal of the matrix for better orthogonality
        	for(int ind = 0; ind< dimension; ind ++)
        	{
            		mastermatrixComplex.data2D[ind][ind] -=1.0;
        	}

		 eigen esolver(1);
		 esolver.GetEigenSolution(mastermatrixComplex, evals, evecsL,  evecsR);
                 
                 // shifting the eigenvalues back
                 cout<<"# Decay constants:\n";
        	for(int ind = 0; ind< dimension; ind ++)
        	{
            	evals.data1D[ind] +=1.0;
            	cout<<evals.data1D[ind]<<"\n";
        	}
                 
        
//        cout<<"evecsR:\n";
//        for(int ind = 0; ind< dimension; ind ++)
//        for(int ins = 0; ins< dimension; ins ++)
//        {
//            cout<<evecsR.data2D[ind][ins];
//            if(ins == dimension-1)
//                cout<<"\n";
//            else
//                cout<<"\t";
//        }
//        cout<<"evecsL:\n";
//        for(int ind = 0; ind< dimension; ind ++)
//            for(int ins = 0; ins< dimension; ins ++)
//            {
//                cout<<evecsL.data2D[ind][ins];
//                if(ins == dimension-1)
//                    cout<<"\n";
//                else
//                    cout<<"\t";
//            }

		// finding largest eigenvalue
		// for proper setup it must be equal to zero
		double max = -1e64;
		int maxindex;
		for(int ind = 0; ind< dimension; ind ++)
		{
			if(max < evals.data1D[ind].real())
			{
				max = evals.data1D[ind].real();
				maxindex = ind;
			}
		}
		if(fabs(max) > constants::smallepsilon)
		{
			cout<<"Error: propagatorM eigenvalue solution failed.\n";
			errcode = 1;
			return retvals;
		}
		else
		{
			// setting strong zero
			//evals.data1D[maxindex] = 0.0;
		}

		pevals = evals.data1D;
		pevecsL = evecsL.data2D;
		pevecsR = evecsL.data2D;

	 }

	 // perform exponentiation and propagation
	//for(int i3=0;i3<dimension; i3++)
	//{
	int& i3 = statefinal;
		retvals = 0.0;
		 for(int i2=0;i2<dimension; i2++)
			 for(int i1=0;i1<dimension; i1++)
				 retvals += pevecsR[i3][i2]*exp(pevals[i2]*time)*pevecsL[i1][i2]*configurationinitialcomplex.data1D[i1];
	//}

	 return retvals;
}


// ***********************
// master equation with memory
// derivs function is impossible to create
// since it is usually time depenent and
// this time dependency is not universal.
// instead it is suggested that 
// the time-depenant kernel is prepared by the user
// at each time step and then the class makes
// a single propagation step
// usage :
// Init, Convolute, Update,  and Exit


// this function initializes the algorithm
//void propagatorM::InitWithMemory(
//                                      int numSt, // number of parameters
//                                      int numM  // number of steps in memory
//                                      )
//{
//    dimension = numSt;
//    dataInt = new int[1];
//    dataInt[0]= numM;
//}
//    
//// this function calculates updates current variables when derivs is given
//void propagatorM::UpdateWithMemory(
//                    complexv** DMM, // current values with history
//                    complexv* der,  // current derivative values
//                    double timS
//)
//{
//    int numM = dataInt[0];
//
//    // updating the density matrix
//    // and shifting memory
//    complexv* tdm = DMM[numM-1];
//    for(int itm=numM-1; itm>0; itm--)
//        DMM[itm] = DMM[itm-1];
//    DMM[0] = tdm;
//
//	// and new result
//    for(int ap = 0; ap<dimension; ap++)
//        DMM[0][ap] = DMM[1][ap] + der[ap]*timS;
//        
//        // done
//}
//// this function calculates derivatives when kernel is given
//// first index of the kernel and of the variable is the history index
//void propagatorM::ConvoluteWithMemory(
//                                           complexv*** kernel,
//                                           complexv** DMM,
//                                           complexv* der
//                                           )
//{
//    int numM = dataInt[0];
//    int& numSt = dimension;
//  
//    for(int indl=0;indl<numSt; indl++)
//    {
//        der[indl]=0.0;
//            for(int indt=0;indt<numM; indt++)
//            for(int indr=0;indr<numSt; indr++)
//                der[indl] += kernel[indt][indl][indr]*DMM[indt][indr];
//    }
//}
//
//void propagatorM::ExitWithMemory()
//{
//    dimension = 0;
//    delete[] dataInt;
//    dataInt = 0;
//}
//
//
//
//
//
//
////this function calculates single propagation step/
////(numerical solution of first order system of differential equations)
//// function
//// int derivs(double, double*, double*, int&)
//// must be created somewhere. 
//// It calculated the right-side derivative values. 
//// Its parameteres:
//// double - cuurent time
//// double* - a set of parameter values (input)
//// double* - a set of derivative values (return)
//// int - a number of variables
//// Function odeStep calculates one single step 
//// or a markovian ODE system
//int propagatorM::odeStep(
//double  time,    // current time befor the next step
//double* data,   // initial values,  on output replaced by new values
//double  step,    // requested step [time dimension]
////int derivs(double, double*, double*, int&), // derivative function
//double* data_bank,   // for temporary data storage: size N*5
//int& N          // number of variables
//)
//{
//	// data initially has initial values
//	// final values are being created
//    
//    //for interpolation
//	double xx[3]={1.0,0.5,0.25};
//	double yy[3];
//	double err;
//	double rezr;
//
//	double step2 = step/2.0;
//	double step4 = step/4.0;
//    
//    //Creating temporary matrices
//	double *Gab, *Gab1, *Gab2, *Gab3, *Gab4, *Gab5, *Gab6, *Gab7;
//	double *GNn;//for derivatives
//    
//	Gab = data;
//	Gab1= &data_bank[0*N];
//	Gab2= &data_bank[1*N];
//	Gab3= Gab1;
//	Gab4= Gab2;
//	Gab5= &data_bank[2*N];
//	Gab6= Gab3;
//	Gab7= &data_bank[3*N];
//	GNn = &data_bank[4*N];
//    
//    //I will do three point interpolation
//    
//    
//    //1-st, 5-th and 7-th substeps
//	RightSideMaster(time,Gab,GNn,N);
//	for(int i=0;i<N;i++)
//	{
//		Gab1[i]=Gab[i]+GNn[i]*step4;
//		Gab5[i]=Gab[i]+GNn[i]*step2;
//		Gab7[i]=Gab[i]+GNn[i]*step;
//	}
//    
//    //2-nd substep
//	RightSideMaster(time+step4,Gab1,GNn,N);
//	for(int i=0;i<N;i++)
//		Gab2[i]=Gab1[i]+GNn[i]*step4;
//    
//    //3-rd substep
//	RightSideMaster(time+2*step4,Gab2,GNn,N);
//	for(int i=0;i<N;i++)
//		Gab3[i]=Gab2[i]+GNn[i]*step4;
//    
//    //4-th substep
//	RightSideMaster(time+3*step4,Gab3,GNn,N);
//	for(int i=0;i<N;i++)
//		Gab4[i]=Gab3[i]+GNn[i]*step4;
//    
//    //6-th substep
//	RightSideMaster(time+step2,Gab5,GNn,N);
//	for(int i=0;i<N;i++)
//		Gab6[i]=Gab5[i]+GNn[i]*step2;
//    
//    
//    //Now all steps are done and I have to interpolate.
//    //Interpolation is done using standart procedure from
//    //Numerical recipies
//    
//    double merror;
//    int retval = 0;
//    toolsInterpolate obj;
//
//	for(int i=0;i<N;i++)
//	{
//		yy[0]=Gab7[i];
//		yy[1]=Gab6[i];
//		yy[2]=Gab4[i];
//        
//        // estimating the measure of error tolerance:
//        merror = max( fabs(yy[0]-yy[1]), fabs(yy[1]-yy[2]) );
//        if(merror == 0.0)
//            rezr = yy[2];
//        else
//        {
//            obj.ratint(xx,yy,3,0.0,rezr,err);
//            if(fabs(err)>merror)
//            {
//                obj.polint(xx,yy,3,0.0,rezr,err);
//                retval = 1;
//            }
//            else if(fabs(err)>merror)
//            {
//                rezr = yy[2];
//                retval = 2;
//            }
//        }
//		Gab[i] = rezr;
//	}
//    
//    return retval;
//}

double propagatorM::max(double a , double b)
{
    return (a>b?a:b);
}
