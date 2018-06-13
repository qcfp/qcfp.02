#pragma once

#include<cmath>


#include<iostream>
#include<fstream>
#include<string>

using namespace std;



#include"../storage/storage.hpp"
#include"../toolsInterpolate/toolsInterpolate.hpp"
#include"../toolsInterpolate2d/toolsInterpolate2d.hpp"
#include"../complexv/complexv.hpp"


// Interpolated linear function based on the dataset:
// main part is stored in the array $values of arbitrary size
// Argument in this part goes from 0 to $xmax and the data is sampled using num_sampl points.
// After $xmax the function (argument is $x)  is zero.

template <class valT> 
class interpolationF2d
{
    //friend class asymptoticLF_complexv;
    
    
  public:


	interpolationF2d();
	interpolationF2d(const interpolationF2d&);
        interpolationF2d(storage<valT>& idat,double ix,double dx,double iy,double dy);
        interpolationF2d(valT** idat,double ix,double dx,int numpx,double iy,double dy,int numpy);
        interpolationF2d(valT* idat,double ix,double dx,int numpx,double iy,double dy,int numpy);
        interpolationF2d(double ix,double dx,int numpx,double iy,double dy,int numpy);
	~interpolationF2d();

	//Returns value of the function for arbitrary argument
	valT Get(double x,double y);

	//operators
	// assignment
	const interpolationF2d& operator = (valT rhs);// this is meaningful 
	//interpolationF& operator = (int rhs);// this is meaningful only when rhs = 0
	//interpolationF& operator = (double rhs);// this is meaningful only when rhs = 0
	const interpolationF2d& operator = (const interpolationF2d& rhs);
        
	const interpolationF2d& operator += (const valT rhs);
	const interpolationF2d& operator -= (const valT rhs);
	const interpolationF2d& operator *= (const valT rhs);
	const interpolationF2d& operator /= (const valT rhs);
	const interpolationF2d& operator += (const interpolationF2d& rhs);
	const interpolationF2d& operator -= (const interpolationF2d& rhs);
	const interpolationF2d& operator *= (const interpolationF2d& rhs);
	const interpolationF2d& operator /= (const interpolationF2d& rhs);


	// multiplication by a scalar
	const interpolationF2d operator * (const valT& rhs) const;
	const interpolationF2d operator / (const valT& rhs) const;
	const interpolationF2d operator + (const valT& rhs) const;
	const interpolationF2d operator - (const valT& rhs) const;

	// linear operations
	const interpolationF2d operator + (const interpolationF2d& rhs) const;
	const interpolationF2d operator - (const interpolationF2d& rhs) const;
	const interpolationF2d operator * (const interpolationF2d& rhs) const;
	const interpolationF2d operator / (const interpolationF2d& rhs) const;


	// comparison
	const bool operator == (const interpolationF2d& rhs);
	const bool operator != (const interpolationF2d& rhs);


        
	int SaveF(string);
	int ReadF(string);
	int SaveF(ofstream*);
	int ReadF(ifstream*);
	int SaveF(ofstream&);
	int ReadF(ifstream&);
        
        
        // accessory functions
        void DirectAssign(interpolationF2d<valT>& irep);
        valT** DirectExchange(valT** irep);
        valT** DirectAccessD();
        void UpdateAxis(double iminx, double istepx,double iminy, double istepy)
        {
            xmin = iminx;
            xstep = istepx;
            xmax = xmin +num_samplx*xstep;
            ymin = iminy;
            ystep = istepy;
            ymax = ymin +num_samply*ystep;
        }
        //Returns value of the function for arbitrary argument
	int GetNx()
        {
            return num_samplx;
        }
	int GetNy()
        {
            return num_samply;
        }
        
        double GetXI()
        {
            return xmin;
        }
        double GetYI()
        {
            return ymin;
        }

        double GetXF()
        {
            return xmax;
        }
        double GetYF()
        {
            return ymax;
        }

  protected:
      	storage<valT>* values;

        int num_samplx;
        int num_samply;

	//step size of argument during initial part
	double xmax;
        double xmin;
        double xstep;
        
	double ymax;
        double ymin;
        double ystep;
        
        
        
        bool CompareType(const interpolationF2d& rhs);
        
        
    private: 
            
        void UpdateMore(double iminx, double istepx,double iminy, double istepy,int numx, int numy)
        {
            num_samplx = numx;
            num_samply = numy;
            xmin = iminx;
            xstep = istepx;
            xmax = xmin +num_samplx*xstep;
            ymin = iminy;
            ystep = istepy;
            ymax = ymin +num_samply*ystep;
        }
        valT GetNegative(double x, double y);


};



//////////////////////////////
// now functions of templates

template<class valT> 
interpolationF2d<valT>::interpolationF2d()
{
	num_samplx = 0;
	num_samply = 0;
	values = new storage<valT>(2);

        xmax = 0.0;
        xmin = 0.0;
        xstep = 0.0;
        ymax = 0.0;
        ymin = 0.0;
        ystep = 0.0;
}

template<class valT> 
interpolationF2d<valT>::interpolationF2d(storage<valT>& idat,double ix,double dx,double iy,double dy)
{
	idat.GetSize(num_samplx,num_samply);
	values = new storage<valT>(2);
	values->Allocate(num_samplx,num_samply);
        for(int ib=0; ib<num_samplx;ib++)
        	values->SetBar(ib,idat.data2D[ib],num_samply);
	
        xmin = ix;
        xstep = dx;
	xmax = ix+xstep*num_samplx;

        ymin = iy;
        ystep = dy;
	ymax = iy+ystep*num_samply;
}
template<class valT> 
interpolationF2d<valT>::interpolationF2d(double ix,double dx,int numpx,double iy,double dy,int numpy)
{
	num_samplx = numpx;
	num_samply = numpy;
	values = new storage<valT>(2);
        values->Allocate(num_samplx,num_samply);
	
        xmin = ix;
        xstep = dx;
	xmax = ix+xstep*numpx;
        ymin = iy;
        ystep = dy;
	ymax = iy+ystep*numpy;
}
template<class valT> 
interpolationF2d<valT>::interpolationF2d(valT** idat,double ix,double dx,int numpx,double iy,double dy,int numpy)
{
	num_samplx = numpx;
	num_samply = numpy;
	values = new storage<valT>(2);
        values->Allocate(num_samplx,num_samply);
        for(int ib=0; ib<num_samplx;ib++)
        	values->SetBar(ib,idat[ib],num_samply);
	
        xmin = ix;
        xstep = dx;
	xmax = ix+xstep*numpx;
        ymin = iy;
        ystep = dy;
	ymax = iy+ystep*numpy;

}
template<class valT> 
interpolationF2d<valT>::interpolationF2d(valT* idat,double ix,double dx,int numpx,double iy,double dy,int numpy)
{
	num_samplx = numpx;
	num_samply = numpy;
	values = new storage<valT>(2);
        values->Allocate(num_samplx,num_samply);
        for(int ib=0; ib<num_samplx;ib++)
        	values->SetBar(ib,idat+ib*num_samply,num_samply);
	
        xmin = ix;
        xstep = dx;
	xmax = ix+xstep*numpx;
        ymin = iy;
        ystep = dy;
	ymax = iy+ystep*numpy;

}
template<class valT> 
interpolationF2d<valT>::interpolationF2d(const interpolationF2d<valT>& rhs)
{
	rhs.values->GetSize(num_samplx,num_samply);
	values = new storage<valT>(2);
	values->Allocate(num_samplx,num_samply);
	
        for(int ib=0; ib<num_samplx;ib++)
        	values->SetBar(ib,rhs.values->data2D[ib],num_samply);
	
        xmin = rhs.xmin;
        xstep = rhs.xstep;
	xmax = rhs.xmax;
        ymin = rhs.ymin;
        ystep = rhs.ystep;
	ymax = rhs.ymax;
}
template<class valT> 
interpolationF2d<valT>::~interpolationF2d()
{
        values->Delete();
        delete values;
        num_samplx = 0;
        num_samply = 0;
        xmax = 0.0;
        xmin = 0.0;
        xstep = 0.0;
        ymax = 0.0;
        ymin = 0.0;
        ystep = 0.0;

}

template<class valT> 
bool interpolationF2d<valT>::CompareType(const interpolationF2d<valT>& rhs)
{

	bool condition= this != &rhs;
        if(condition)
        {
                condition = condition && this->xmax==rhs.xmax;
                condition = condition &&  this->xmin==rhs.xmin;
                condition = condition &&  this->xstep==rhs.xstep;
                condition = condition &&  this->num_samplx == rhs.num_samplx;
                condition = condition && this->ymax==rhs.ymax;
                condition = condition &&  this->ymin==rhs.ymin;
                condition = condition &&  this->ystep==rhs.ystep;
                condition = condition &&  this->num_samply == rhs.num_samply;
        }
        else return true;
        
        return condition;
}



template<class valT> 
int interpolationF2d<valT>::SaveF(ofstream& fstr)
{
    SaveF(&fstr);
}
template<class valT> 
int interpolationF2d<valT>::SaveF(ofstream* fstr)
{
    char message[100] = "Interpolation function 2d, v. 2.0\n";

        
	fstr->write(message, 100*sizeof(char));
	fstr->write( (char*) &num_samplx,sizeof(int));
	fstr->write( (char*) &xmin,sizeof(double));
	fstr->write( (char*) &xmax,sizeof(double));
	fstr->write( (char*) &xstep,   sizeof(double));
	fstr->write( (char*) &num_samply,sizeof(int));
	fstr->write( (char*) &ymin,sizeof(double));
	fstr->write( (char*) &ymax,sizeof(double));
	fstr->write( (char*) &ystep,   sizeof(double));
        
	if(num_samplx>0 && num_samply>0) 
                for(int ib=0; ib<num_samplx;ib++)
		fstr->write( (char*) values->data2D[ib],num_samply*sizeof(valT));
}


template<class valT> 
int interpolationF2d<valT>::SaveF(string str)
{
    
    ofstream* fstr = new ofstream(str.c_str());

    SaveF(fstr);      
    
    fstr->close();
        
    delete fstr;
}
template<class valT> 
int interpolationF2d<valT>::ReadF(ifstream& fstr)
{
       ReadF(&fstr);
}
template<class valT> 
int interpolationF2d<valT>::ReadF(string str)
{
       ifstream* fstr = new ifstream(str.c_str());
       ReadF(fstr);
       fstr->close();
       delete fstr;
}
template<class valT> 
int interpolationF2d<valT>::ReadF(ifstream* fstr)
{

	char message[100];

	fstr->read(message, 100*sizeof(char));
        
	cout<<"Reading: "<<string(message);

        fstr->read( (char*) &num_samplx,sizeof(int));
        fstr->read( (char*) &xmin,sizeof(double));
        fstr->read( (char*) &xmax,sizeof(double));
        fstr->read( (char*) &xstep,   sizeof(double));
        fstr->read( (char*) &num_samply,sizeof(int));
        fstr->read( (char*) &ymin,sizeof(double));
        fstr->read( (char*) &ymax,sizeof(double));
        fstr->read( (char*) &ystep,sizeof(double));

	if(num_samplx>0 && num_samply>0) 
	{
                values->Delete();
                values->Allocate(num_samplx,num_samply);
                for(int ib=0; ib<num_samplx;ib++)
                fstr->read( (char*) values->data2D[ib],num_samply*sizeof(valT));
	}
	else
	{
                cout<<"Trying to read an empty function. Aborting.\n";
	}
}


template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator = (const interpolationF2d<valT>& rhs)
{
	if (this != &rhs)
	{
                xmax=rhs.xmax;
                xmin=rhs.xmin;
                xstep=rhs.xstep;
                ymax=rhs.ymax;
                ymin=rhs.ymin;
                ystep=rhs.ystep;

                if(num_samplx != rhs.num_samplx || num_samply != rhs.num_samply )
                {
                    values->Delete();
                    num_samplx = rhs.num_samplx;
                    num_samply = rhs.num_samply;
                    values->Allocate(num_samplx,num_samply);
                }
                for(int ib=0; ib<num_samplx;ib++)
                values->SetBar(ib,rhs.values->data2D[ib],num_samply);

	}
	return *this;    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator = (const valT rhs)
{ 
	{
                xmax=1;
                xmin=0;
                xstep=1;
                ymax=1;
                ymin=0;
                ystep=1;

                    values->Delete();
                    values->Allocate(1,1);
                    num_samplx = 1;
                    num_samply = 1;
        
		// getting lhs data object
		values->data2D[0][0] = rhs;

	}
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator += (const valT rhs)
{ 
	{
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]+=rhs;
        
	}
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator -= (const valT rhs)
{ 
	{        
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]-=rhs;
	}
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator *= (const valT rhs)
{ 
	{
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]*=rhs;
        
	}
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator /= (const valT rhs)
{ 
	{        
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]/=rhs;
	}
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator += (const interpolationF2d<valT>& rhs)
{
        if(CompareType(rhs))
        {        
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]+=rhs.values->data2D[indx][indy];
        }
        else
        {
             cout<<"Error: cannot append object (file interpolationF.hpp)\n";
        }
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator -= (const interpolationF2d<valT>& rhs)
{
        if(CompareType(rhs))
        {        
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]-=rhs.values->data2D[indx][indy];
        }
        else
        {
             cout<<"Error: cannot append object (file interpolationF.hpp)\n";
        }
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator *= (const interpolationF2d<valT>& rhs)
{
        if(CompareType(rhs))
        {        
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]*=rhs.values->data2D[indx][indy];
        }
        else
        {
             cout<<"Error: cannot append object (file interpolationF.hpp)\n";
        }
	return *this;    // Return ref for multiple assignment
}
template<class valT> 
const interpolationF2d<valT>& interpolationF2d<valT>::operator /= (const interpolationF2d<valT>& rhs)
{
        if(CompareType(rhs))
        {        
                for(int indx=0; indx<num_samplx; indx++)
                 for(int indy=0; indy<num_samply; indy++)
       		this->values->data2D[indx][indy]/=rhs.values->data2D[indx][indy];
        }
        else
        {
             cout<<"Error: cannot append object (file interpolationF.hpp)\n";
        }
	return *this;    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator + (const interpolationF2d<valT>& rhs)
 const
{
	return interpolationF2d<valT>(*this += rhs);    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator - (const interpolationF2d<valT>& rhs)
 const
{
	return interpolationF2d<valT>(*this -= rhs);    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator * (const interpolationF2d<valT>& rhs)
 const
{
	return interpolationF2d<valT>(*this *= rhs);    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator / (const interpolationF2d<valT>& rhs)
 const
{
	return interpolationF2d<valT>(*this /= rhs);    // Return ref for multiple assignment
}


template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator + (const valT& rhs)
 const
{
	return interpolationF2d<valT>(*this += rhs);    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator - (const valT& rhs)
 const
{
	return interpolationF2d<valT>(*this -= rhs);    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator * (const valT& rhs)
 const
{
	return interpolationF2d<valT>(*this *= rhs);    // Return ref for multiple assignment
}

template<class valT> 
const interpolationF2d<valT> interpolationF2d<valT>::operator / (const valT& rhs)
 const
{
	return interpolationF2d<valT>(*this /= rhs);    // Return ref for multiple assignment
}


template<class valT> 
const bool interpolationF2d<valT>::operator == (const interpolationF2d<valT>& rhs)
{
    if(this == &rhs)
        return true;
    
    bool result = CompareType(rhs);
    
    if(result)
    {

                for(int indx=0; indx<this->num_samplx; indx ++)
                for(int indy=0; indy<this->num_samply; indy ++)
                        result = result && (this->values->data2D[indx][indy]==rhs.values->data2D[indx][indy]);
    }
        return result;
}

template<class valT> 
const bool interpolationF2d<valT>::operator != (const interpolationF2d<valT>& rhs)
{
	return !(operator ==(rhs));    // Return ref for multiple assignment
}


template<class valT> 
void interpolationF2d<valT>::DirectAssign(interpolationF2d<valT>& irep)
{
    if(num_samplx<irep.num_samplx || num_samply<irep.num_samply )
    {
        cout<<"Error: in void interpolationF2d<valT>::DirectReplace(interpolationF2d<valT>& irep)\n";
        return;
    }
    else
    {
    for(int indx = 0; indx< irep.num_samplx; indx++ )
    for(int indy = 0; indy< irep.num_samply; indy++ )
    {
        values->data2D[indx][indy] = irep.values->data2D[indx][indy];
    }
    }
}
template<class valT> 
valT** interpolationF2d<valT>::DirectExchange(valT** irep)
{
    valT** tmp = values->data2D;
    values->data2D = irep;
    return tmp;
}


template<class valT> 
valT** interpolationF2d<valT>::DirectAccessD()
{
    return values->data2D;
}

template<class valT> 
valT interpolationF2d<valT>::GetNegative(double x, double y)
{
    if(xstep<0 && ystep<0)
    {
        double distx = xmin - x;
        double disty = ymin - y;
        interpolationF2d<valT> tfunct;
        valT** tptr = tfunct.DirectExchange(values->data2D);
        tfunct.UpdateMore(xmin, -xstep,ymin, -ystep,num_samplx, num_samply);
        valT retval = tfunct.Get(xmin+distx,ymin+disty);
        values->data2D = tfunct.DirectExchange(tptr);
        return retval;
    }
    else if(xstep<0 && ystep>0)
    {
        double distx = xmin - x;
        //double disty = ymin - y;
        interpolationF2d<valT> tfunct;
        valT** tptr = tfunct.DirectExchange(values->data2D);
        tfunct.UpdateMore(xmin, -xstep,ymin, ystep,num_samplx, num_samply);
        valT retval = tfunct.Get(xmin+distx,y);
        values->data2D = tfunct.DirectExchange(tptr);
        return retval;
    }
    if(xstep>0 && ystep<0)
    {
        //double distx = xmin - x;
        double disty = ymin - y;
        interpolationF2d<valT> tfunct;
        valT** tptr = tfunct.DirectExchange(values->data2D);
        tfunct.UpdateMore(xmin, xstep,ymin, -ystep,num_samplx, num_samply);
        valT retval = tfunct.Get(x,ymin+disty);
        values->data2D = tfunct.DirectExchange(tptr);
        return retval;
    }
    return (valT)0.0;
 
    cout<<"Error in interpolationF2d<valT>::GetNegative(double x, double y)\n";
}