#pragma once


#include"../storage/storage.hpp"
#include"../complexv/complexv.hpp"

class propagatorM
{
public:
    
    // simple master equation
    // uses eigenvalues and eigenvectors for propagation
    // rate matrix must be real-valued

	propagatorM();
	~propagatorM();
    // usage 
    // Init, Propagate and Exit
    int InitMaster(storage<double>& matrix);
    int InitMaster(storage<complexv>& matrix);
//    int InitMaster(int num, double** matrix);
//    int PropagateMaster(double* y,int numt,double dt);
    int PropagateMaster(storage<double>& y,int numt,double dt);
    int ExitMaster();

    // returns exponential  starting from initial conditions
    // y
    storage<double> Get(storage<double>& y,double& time);
    storage<complexv> Get(storage<complexv>& y,double& time);
    
    // returns exponential for at specific state,
    // starting from stores initials
    double GetAt(double& time);
    complexv GetAtC(double& time);
    void prepareGetAt(storage<double>& initials,int retindex);
    void prepareGetAt(storage<complexv>& initials,int retindex);
    
//    // this function initializes the algorithm
//    void InitWithMemory(int numSt, // number of parameters
//                        int numM  // number of steps in memory
//                        );
    
//    // this function calculates updates current variables when derivs is given
//    void UpdateWithMemory(complexv** DMM, // current values with history
//                          complexv* der,  // current derivative values
//                          double dt // time step
//                          );
    
//    // this function calculates derivatives when kernel is given
//    void ConvoluteWithMemory(complexv*** kernel,
//                             complexv** DMM,
//                             complexv* der
//                             );
//    void ExitWithMemory();
    
    
    
    // universal function calculates single propagation step for ODE
    // so it is OK for the Markovian differential equations
    //(numerical solution of first order system of differential equations)
    // function
    // int derivs(double, double*, double*, int&)
    // must be created somewhere. 
    // It calculated the right-side derivative values. 
    // Its parameteres:
    // double - cuurent time
    // double* - a set of parameter values (input)
    // double* - a set of derivative values (return)
    // int - a number of variables
    // Function odeStep calculates one single step 
    // or a markovian ODE system
//    int odeStep(
//                double time,    // current time befor the next step
//                double* data,   // initial values,  on output replaced by new values
//                double step,    // requested step [time dimension]
//                //int derivs(double, double*, double*, int&),
//                             // derivative function
//                double* data_bank,   // for temporary data storage: size N*5
//                int& N          // number of variables
//                );
//

public:
    complexv* pevals;
    complexv** pevecsR;
    complexv** pevecsL;
    storage<complexv> evals;
    storage<complexv> evecsR;
    storage<complexv> evecsL;

    
private:
    
    double max(double a , double b);

//    int RightSideMaster(double t, double* y, double* dy, int& num);
    
    int dimension;
    int* dataInt;
    //storage<double> dataset1D;
    storage<double> mastermatrix;
    storage<complexv> mastermatrixComplex;
    
    storage<double> configurationinitial;
    storage<complexv> configurationinitialcomplex;
    int statefinal;



	int errcode;
	// nonzero value means that propagation is failed

};

