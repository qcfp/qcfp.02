#pragma once


#include"../storage/storage.hpp"
#include"../interpolationF/interpolationF.hpp"



template<class valT, class derivsT>
class propagatorODEa
{
    // interface:
public:
    propagatorODEa(derivsT* isystemObj, int numDF);
    
    int propagateODEstep(
                          double time,    // current time before the next step
                          valT* iodata,   // initial values,  on output replaced by new values
                          double step     // requested step [time dimension]
    );
    storage<interpolationF<valT> > propagateODE(
                                     storage<double>& time,  // all times
                                     storage<valT>& iodata   // initial values,
    );

	int GetNumDF()
	{
		return numDF;
	}


private:
    storage<valT> data_bank; // for temporary data storage: size Num*5
    derivsT* systemObj;
    int numDF;
    
    valT *Gab1, *Gab2, *Gab3, *Gab4, *Gab5, *Gab6, *Gab7;
    valT *GNn;//for derivatives
    
    

};


/////////////////////////////
// implementations
/////////////////////////////

template<class valT, class derivsT>
propagatorODEa< valT,  derivsT>::propagatorODEa(derivsT* isystemObj, int inumDF)
{
    systemObj = isystemObj;
    numDF = inumDF;
    
    data_bank.SetDimension(1);
    data_bank.Allocate(5*numDF);
    
    Gab1= data_bank.data1D;
    Gab2= data_bank.data1D+numDF;
    Gab3= Gab1;
    Gab4= Gab2;
    Gab5= data_bank.data1D+2*numDF;
    Gab6= Gab3;
    Gab7= data_bank.data1D+3*numDF;
    GNn = data_bank.data1D+4*numDF;

}

// format without storage class
template<class valT, class derivsT>
int propagatorODEa< valT,  derivsT>::propagateODEstep(double time, valT* iodata, double step)
{
    // data initially has initial values
    // final values are being created
    int& N = numDF;
    valT* state;
    valT* deriv;
    
    double step2 = step/2.0;
    double step4 = step/4.0;
    
    //I will do three point interpolation
    
    //1-st, 5-th and 7-th substeps
    state = iodata;
    systemObj->operator()(time,state,GNn);
    for(int i=0;i<N;i++)
    {
        Gab1[i]=state[i]+GNn[i]*step4;
        Gab5[i]=state[i]+GNn[i]*step2;
        Gab7[i]=state[i]+GNn[i]*step;
    }
    
    //2-nd substep
    state = Gab1;
    systemObj->operator()(time+step4,state,GNn);
    for(int i=0;i<N;i++)
        Gab2[i]=state[i]+GNn[i]*step4;
    
    //3-rd substep
    state = Gab2;
    systemObj->operator()(time+step2,state,GNn);
    for(int i=0;i<N;i++)
        Gab3[i]=state[i]+GNn[i]*step4;
    
    //4-th substep
    state = Gab3;
    systemObj->operator()(time+step2+step4,state,GNn);
    for(int i=0;i<N;i++)
        Gab4[i]=state[i]+GNn[i]*step4;
    
    //6-th substep
    state = Gab5;
    systemObj->operator()(time+step2,state,GNn);
    for(int i=0;i<N;i++)
        Gab6[i]=state[i]+GNn[i]*step2;
    
    
    //Now all steps are done and I have to interpolate.
    
    // making the approaching differences
    
    for(int i=0;i<N;i++)
    {
        valT dif2 = (Gab7[i]-Gab6[i])/step2;
        valT dif1 = (Gab6[i]-Gab4[i])/step4;
        
        valT dif0 = 2.0*dif1 - dif2;
        
        // the return value
        iodata[i] = Gab4[i] - dif0*step4;
    }
    
    return 0;
    
}

template<class valT, class derivsT>
storage<interpolationF<valT> > propagatorODEa< valT,  derivsT>::propagateODE(storage<double>& time, storage<valT>& iodata)
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
