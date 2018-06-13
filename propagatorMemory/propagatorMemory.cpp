
#include "propagatorMemory.hpp"
#include"../toolsInterpolate/toolsInterpolate.hpp"
#include "../storage/storage.hpp"
#include "../constants/constants.hpp"
#include "../eigen/eigen.hpp"
#include<cmath>


propagatorMemory::propagatorMemory()
:propagatorM()
{
    
    omegas_memory.SetDimension(1);
    timesinternal.SetDimension(1);
    
    kernel= 0;
    DMM=0;
    
    omegas_memory1.SetDimension(1);
    kernel1= 0;
    DMM1=0;
    
    omegas_memory2.SetDimension(1);
    kernel2= 0;
    DMM2=0;
    
    manifold0End=1e99;
    manifold1End=1e99;
    manifold2End=1e99;
    
}

propagatorMemory::~propagatorMemory()
{
    
}

// ***********************
// master equation with memory
// 
// the translationary invariant time-depenant kernel is prepared by the user
// the class makes
// propagation
// usage :
// Init, Convolute, Update,  and Exit


// this function can change the kernel on a fly
    void propagatorMemory::Initialize(storage<complexv>& ikernel, 
                    storage<double>& iomegas_memory, 
                    storage<complexv>& iDMM, 
                    storage<double>& itimesinternal)
    {
        kernel = &ikernel;
        omegas_memory = iomegas_memory;
        DMM = &iDMM;
        timesinternal = itimesinternal;
        
        int timeN = timesinternal.GetSize();
        manifold0End = timesinternal.data1D[timeN-1]-constants::smallepsilon;
    }
    
    void propagatorMemory::Initialize(propagatorMemory& prop)
    {
        kernel = prop.kernel;
        omegas_memory = prop.omegas_memory;
        DMM = prop.DMM;
        timesinternal = prop.timesinternal;
        
        
        manifold0End = prop.manifold0End;
    }
    
    void propagatorMemory::InitializeHistory(propagatorMemory& prop)
    {
        if(prop.kernel != 0){
            kernel1 = prop.kernel;
            omegas_memory1 = prop.omegas_memory;
            DMM1 = prop.DMM;
        //timesinternal = prop.timesinternal;
            manifold1End = prop.manifold0End;
        }
        

        if(prop.kernel1 != 0){
            kernel2 = prop.kernel1;
            omegas_memory2 = prop.omegas_memory1;   
            DMM2 = prop.DMM1;        
        
            manifold2End = prop.manifold1End;
        }
    }


// this function  updates current variables when derivs is given
void propagatorMemory::Update(
    storage<complexv>& der,  // current derivative values
    double dtimes
)
{
    // updating the density matrix
    // and shifting memory
    double timestep = timesinternal.data1D[0]-timesinternal.data1D[1];
    int timeN = timesinternal.GetSize();
    double timeini = timesinternal.data1D[0];
    int dimension = der.GetSize();
    
    complexv* tdm = DMM->data2D[timeN-1];
    for(int itm=timeN-1; itm>0; itm--)
        DMM->data2D[itm] = DMM->data2D[itm-1];
    DMM->data2D[0] = tdm;
    
    // and new result
    for(int ap = 0; ap<dimension; ap++)
    {
        DMM->data2D[0][ap] = DMM->data2D[1][ap] + der.data1D[ap]*dtimes;
    }
    
    // updating times
    timesinternal.FillLinear(timeini+dtimes,-timestep,timeN);
    // done
}

// this function calculates derivatives when kernel is given
// first index of the kernel and of the variable is the history index
void propagatorMemory::Convolute(
    storage<complexv>& der  // current derivative values
)
{
    double timestep = timesinternal.data1D[0]-timesinternal.data1D[1];
    int timeN = timesinternal.GetSize();
    int dimension = der.GetSize();
    int dimension1 = der.GetSize();
    int dimension2 = der.GetSize();
    //std::cout<<"Error: this function cannot be called at this time because of dimension1 and dimension2 incorrect setting\n";
    
    // integration 
    for(int indl=0;indl<dimension; indl++)
    {
        double omegalxtimel = omegas_memory.data1D[indl]*timesinternal.data1D[0];
        complexv& derdata1Dindl = der.data1D[indl];
        derdata1Dindl=0.0;
        
        int indt1=0;
        int indt2=0;
        int indt0=0;
        
        for(int indt=0;indt<timeN; indt++)
        {
            // notice: time goes in reverse direction
            double time = timesinternal.data1D[indt];
            
            
            // convoluting in the same manifold
            if(time > manifold0End)// here one must include the small parameter
            {
                // convolution with itself
                for(int indr=0;indr<dimension; indr++)
                {
                    double omega = omegas_memory.data1D[indr];
                    derdata1Dindl += kernel->data3D[indt][indl][indr]*exp(coni*(omegalxtimel-omega*time))*DMM->data2D[indt0][indr]*timestep;
                }
                indt0++;
            }
            else if(time > manifold1End)
            {
                // convolution with previous manifold
                for(int indr=0;indr<dimension1; indr++)
                {
                    double omega = omegas_memory1.data1D[indr];
                    derdata1Dindl += kernel1->data3D[indt][indl][indr]*exp(coni*(omegalxtimel-omega*time))*DMM1->data2D[indt1][indr]*timestep;
                }
                indt1++;
            }
            else if(time > manifold2End)
            {
                // convolution with even earlier manifold
                for(int indr=0;indr<dimension2; indr++)
                {
                    double omega = omegas_memory2.data1D[indr];
                    derdata1Dindl += kernel2->data3D[indt][indl][indr]*exp(coni*(omegalxtimel-omega*time))*DMM2->data2D[indt2][indr]*timestep;
                }
                
                indt2++;
            }
            
            
            
        }
    }
}

// propagates DMM ]ith respect to itself and history DMM1 and DMM2
storage<complexv> propagatorMemory::Propagate(storage<double>& times)
// historyin DMM: 
// 0 latest time
// 1, 2, ... timeN - deeper and deeper into history

// time 0 corresponds to the currenttime

// time in ret as requested by times

{
    
    // external propagation parameters
    double tdext = times.data1D[1]-times.data1D[0];
    double tiext = times.data1D[0];
    int tvalsext = times.GetSize();
    
    // internal propagation parameters
    double tdint = timesinternal.data1D[0]-timesinternal.data1D[1];
    int numT,dimension;
    DMM->GetSize(numT,dimension);
    
    
    storage<complexv> ret(2); 
    ret.Allocate(numT,dimension);
    
    storage<complexv> derivs(1);
    derivs.Allocate(dimension);
    
    //loops over external times
    for(int itext=0;itext<tvalsext;itext++)
    {
        double timeext = times.data1D[itext];
        
        
        // propagation with small internal timesteps;
        int numloops = (int)((timeext-timesinternal.data1D[0])/tdint); // this should return floor value
        double reminder = timeext-timesinternal.data1D[0]-numloops*tdint;
        
        for(int it=0;it<numloops;it++)
        {
            Convolute( derivs);
            Update(derivs,tdint);
        }
        
        
        // doing reminder assignment:
        Convolute( derivs);
        // and new result
        for(int ap = 0; ap<dimension; ap++)
        {
            double omega = omegas_memory.data1D[ap];
            ret.data2D[itext][ap] = exp(cnni*omega*timeext)*(DMM->data2D[0][ap] + derivs.data1D[ap]*reminder);
        }
    }
    
    return ret;
}


//accessory functions:
storage<complexv> propagatorMemory::Convert2DTo1D(storage<complexv> input)
{
    storage<complexv> result(1);
    int numl,numr;
    input.GetSize(numl,numr);
    result.Allocate(numl*numr);
    for(int idl = 0; idl<numl; idl++)
        for(int idr = 0; idr<numr; idr++)
        {
            result.data1D[numr*idl+idr] = input.data2D[idl][idr];
        }
        return result;
    
}

storage<complexv> propagatorMemory::Convert1DTo2D(storage<complexv> input,int numl,int numr)
{
    storage<complexv> result(2);
    result.Allocate(numl,numr);
    for(int idl = 0; idl<numl; idl++)
        for(int idr = 0; idr<numr; idr++)
        {
            result.data2D[idl][idr]=input.data1D[numr*idl+idr] ;
        }
        return result;
}
storage<complexv> propagatorMemory::Convert2DTo3D(storage<complexv> input,int numl,int numr)
{
    int num2,num1;
    input.GetSize(num2,num1);
    storage<complexv> result(3);
    result.Allocate(num2,numl,numr);
    for(int it = 0; it<num2; it++)
        for(int idl = 0; idl<numl; idl++)
            for(int idr = 0; idr<numr; idr++)
            {
                result.data3D[it][idl][idr]=input.data2D[it][numr*idl+idr] ;
            }
            return result;
}

storage<complexv> propagatorMemory::Convert4DTo2D(storage<complexv> input)
{
    storage<complexv> result(2);
    int num4,num3,num2,num1;
    input.GetSize(num4,num3,num2,num1);
    result.Allocate(num4*num3,num2*num1);
    for(int idl2 = 0; idl2<num4; idl2++)
        for(int idr2 = 0; idr2<num3; idr2++)
            for(int idl = 0; idl<num2; idl++)
                for(int idr = 0; idr<num1; idr++)
                {
                    result.data2D[num3*idl2+idr2][num1*idl+idr] = input.data4D[idl2][idr2][idl][idr];
                }
                return result;
}
storage<complexv> propagatorMemory::Convert2DTo4D(storage<complexv> input,int num4,int num3,int num2,int num1)
{
    storage<complexv> result(4);
    result.Allocate(num4,num3,num2,num1);
    for(int idl2 = 0; idl2<num4; idl2++)
        for(int idr2 = 0; idr2<num3; idr2++)
            for(int idl = 0; idl<num2; idl++)
                for(int idr = 0; idr<num1; idr++)
                {
                    result.data4D[idl2][idr2][idl][idr]=input.data2D[num3*idl2+idr2][num1*idl+idr];
                }
                return result;
}

storage<complexv> propagatorMemory::Convert5DTo3D(storage<complexv> input)
{
    storage<complexv> result(3);
    int num5,num4,num3,num2,num1;
    input.GetSize(num5,num4,num3,num2,num1);
    result.Allocate(num5, num4*num3,num2*num1);
    for(int it = 0; it<num5; it++)
        for(int idl2 = 0; idl2<num4; idl2++)
            for(int idr2 = 0; idr2<num3; idr2++)
                for(int idl = 0; idl<num2; idl++)
                    for(int idr = 0; idr<num1; idr++)
                    {
                        result.data3D[it][num3*idl2+idr2][num1*idl+idr] = input.data5D[it][idl2][idr2][idl][idr];
                    }
                    return result;
}

storage<complexv> propagatorMemory::Convert3DTo5D(storage<complexv> input,int num4,int num3,int num2,int num1)
{
    storage<complexv> result(5);
    
    int num5,tnum2,tnum1;
    input.GetSize(num5,tnum2,tnum1);
    
    
    result.Allocate(num5,num4,num3,num2,num1);
    for(int it = 0; it<num5; it++)
        for(int idl2 = 0; idl2<num4; idl2++)
            for(int idr2 = 0; idr2<num3; idr2++)
                for(int idl = 0; idl<num2; idl++)
                    for(int idr = 0; idr<num1; idr++)
                    {
                        result.data5D[it][idl2][idr2][idl][idr]=input.data3D[it][num3*idl2+idr2][num1*idl+idr];
                    }
                    return result;
}



