/* 
 * File:   propagatorODErk.hpp
 * Author: dariusa
 *
 * Created on October 27, 2015, 5:34 PM
 */

#ifndef PROPAGATORODERK_HPP
#define	PROPAGATORODERK_HPP

// runge-kutta propagator




#pragma once


#include"../storage/storage.hpp"
#include"../interpolationF/interpolationF.hpp"



template<class valT, class derivsT>
class propagatorODErk
{
    // interface:
public:
    propagatorODErk(derivsT* isystemObj, int numDF);
    
    int propagateODEstep(
                          double time,    // current time before the next step
                          valT* iodata,   // initial values,  on output replaced by new values
                          double step     // requested step [time dimension]
    );
    storage<interpolationF<valT> > propagateODE(
                                     storage<double>& time,  // all times
                                     storage<valT>& iodata   // initial values,
    );



private:
    storage<valT> data_bank; // for temporary data storage: size Num*5
    derivsT* systemObj; // stores the equations object
    int numDF;
    
    valT *yk,*k1,*k2,*k3,*k4;
    
    

};


/////////////////////////////
// implementations
/////////////////////////////

template<class valT, class derivsT>
propagatorODErk< valT,  derivsT>::propagatorODErk(derivsT* isystemObj, int inumDF)
{
    systemObj = isystemObj;
    numDF = inumDF;
    
    data_bank.SetDimension(1);
    data_bank.Allocate(5*numDF);
    
    yk= data_bank.data1D;
    k1= data_bank.data1D+numDF;
    k2= data_bank.data1D+2*numDF;
    k3= data_bank.data1D+3*numDF;
    k4= data_bank.data1D+4*numDF;

}

// format without storage class
template<class valT, class derivsT>
int propagatorODErk< valT,  derivsT>::propagateODEstep(double time, valT* iodata, double step)
{
    // data initially has initial values
    // final values are being created
    int& N = numDF;
    valT* state;
    valT* deriv;
    
    double step2 = step/2.0;
    double step6 = step/6.0;
    
 
    
    
    
    //1-st, substep
    systemObj->operator()(time,iodata,k1);
    for(int i=0;i<N;i++)
        yk[i]=iodata[i]+k1[i]*step2;
    
    //2-nd substep
    systemObj->operator()(time+step2,yk,k2);
    for(int i=0;i<N;i++)
        yk[i]=iodata[i]+k2[i]*step2;
    
    //3-rd substep
    systemObj->operator()(time+step2,yk,k3);
    for(int i=0;i<N;i++)
        yk[i]=iodata[i]+k3[i]*step;
    
    //4-th substep
    systemObj->operator()(time+step,yk,k4);

    for(int i=0;i<N;i++)
        iodata[i]=iodata[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])*step6;
    

    return 0;
    
}

template<class valT, class derivsT>
storage<interpolationF<valT> > propagatorODErk< valT,  derivsT>::propagateODE(storage<double>& time, storage<valT>& iodata)
{
    
    int numT = time.GetSize();
    double tstep;
    if(numT>1)
        tstep = time.data1D[1]-time.data1D[0];
    else tstep = 0;

    storage<valT> dataset(2);
    dataset.Allocate(numDF,numT);
    
    for(int is = 0; is< numDF; is++)
        dataset.data2D[is][0] =iodata.data1D[is];
 
    for(int it=0; it<numT-1; it++)
    {
        double tini= time.data1D[it];
        propagateODEstep(tini, iodata.data1D, tstep);
        
        for(int is = 0; is< numDF; is++)
            dataset.data2D[is][it+1] =iodata.data1D[is];

    }
    
    storage<interpolationF<valT> > retval(1);
    retval.Allocate(numDF);
    
    for(int is = 0; is< numDF; is++)
        retval.data1D[is] = interpolationF<valT>(dataset.data2D[is],time.data1D[0],tstep,numT);

    
    return retval;
}



#endif	/* PROPAGATORODERK_HPP */

