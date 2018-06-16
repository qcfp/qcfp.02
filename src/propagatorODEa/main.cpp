// this is the example program to solve the van der pol equation



#include"../propagatorODEa/propagatorODEa.hpp"
#include <iomanip>      // std::setprecision



class vanderpol
{
public:
    // parameters
    double mu;
    double omega;
    
public:
    
    vanderpol(double imu, double iomega)
    {
        // two parameters
        mu = imu;
        omega = iomega;
    }
    
    
    void operator ()(double itime, double* state_now,  double* derivatives)
    // returns an array of derivative values based on the ODE
    {
        double& x = state_now[0];
        double& y = state_now[1];
        derivatives[0]= y;
        derivatives[1]= -omega*omega*x +mu*(1-x*x)*y;
    }
    
};


int main()
{
    
    vanderpol objvdp(1, 2*3.1416);
    
    propagatorODEa<double,vanderpol> objprop(&objvdp, 2);
    
    
    // internal current state
    storage<double> state_current;
    // two degrees of freedom
    state_current.SetDimension(1);
    state_current.Allocate(2);
    
    state_current.data1D[0] = 1;
    state_current.data1D[1] = 0.0;

    
    storage<double> time;
    double step = 0.1;
    int nump=100;
    time.FillLinear(0, step, nump);

    storage<interpolationF<double> > solution = objprop.propagateODE(time, state_current);
    
    // making out dataset for plotting
    for(int it = -2; it<nump*10+2; it++)
    {
        double ltime =(step*it)/10;
        cout<<ltime<<"\t"<<solution.data1D[0].Get(ltime)<<"\t"<<solution.data1D[1].Get(ltime)<<"\n";
    }

    
}