// Asymptotic linear function in numeric format
// the part from xmin to xmax-step is stored explicitly in numerical format
// using a step size of step.
// the number of points is as many as int allows
// values between given points are interpolated using Splines
// outside the interval the function is assumed to be linear function

// The numerical part is stored in the array $values
// Argument in this part goes from 0 to $xmax and the data is sampled using num_sampl points.

// for arguments larger or equal than xmax the function (argument is x)  is linear: a*(x-xmax)+b.
// for arguments lower or equal than xmin the function is linear: c*(x-xmin)+d


/////////////////
// this file defines

// template for an arbitrary data types
// asymptoticLF

// specific float type (typedef)
// asymptoticLF_float


// specific double type (typedef)
// asymptoticLF_double

// specific complexv type (different than typedef: see below)
// asymptoticLF_complexv


#pragma once


//using namespace std;


#include<iostream>
#include<fstream>
#include<string>
#include"../storage/storage.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../interpolationF/interpolationF.hpp"
#include"../complexv/complexv.hpp"




template <class valT>
	class asymptoticLF:
public interpolationF<valT>
{
  public:

      

	//asymptotic linear coefficients
	//starting from xmax till infinity the function is valT(a*(argument-xmax))+b
	valT a;
    // a is a derivative of the function at xmax
	valT b;
    //b=function(xmax)

	valT c;
    // c is the derivative at xmin
	valT d;
	//starting from -infinity till xmin the function is valT(c*(argument-xmin))+d
    // d is the value at xmin
    // d=val[0]
        
        valT deltaamplitude;
        // at zero time there is additional dirac delta function
        // this is a factor in front of the dirac delta.
        
    // special flag about symmetry:
    int causality;
    // 0 means ignored
    // 2 means  f(-t) = f(t)
    // 3 means  f(-t) = -f(t)

    
    asymptoticLF();
    // simple empty constructor
    
	asymptoticLF(const asymptoticLF<valT>& idat);
    // a copy is made
    
	asymptoticLF(storage<valT>& idat,double ixmin,double istep);
    // a function is created from a dataset. Derivatives are calculated from the dataset
    // xmin and step are given
    // number of points is provided in the storage container
    
	asymptoticLF(storage<valT>& idat,double ixmin,double istep, valT iderivI, valT ivalF, valT iderivF);
    // a function is created from a dataset. Derivatives are provided:
    // iderivI - derivative at xmin
    // iderivF - derivative at xmax
    // ivalF - function value at xmax

        
    asymptoticLF(valT* idat, double ixmin,double istep,int inump);
    // a function is created from a linear array. Derivatives are calculated from the dataset
    // xmin and step are given
    // number of points is given

	~asymptoticLF();
    // standard destructor

	//Returns value of the function for an arbitrary argument
        // this does not include the dirac delta at zero
        // next function returns the factor of delta
	valT Get(double x);

	valT GetDiracAmplitude()
        {
            if(causality == 3)
                return 0;
            else return deltaamplitude;
        };
	void SetDiracAmplitude(valT ival)
        {
            if(causality != 3)
                deltaamplitude = ival;
        };

        
        // returns the pure zero element
    valT Get(string& x)
    {
        if(interpolationF<valT>::num_sampl != 0)
            return  interpolationF<valT>::values->data1D[0];
        else return 0.0;
    }

    // standard operators on the function
	const asymptoticLF& operator = (const valT rhs);// this is meaningful 
	const asymptoticLF& operator = (const asymptoticLF& rhs);
	const asymptoticLF& operator = ( interpolationF<valT>& rhs);
        
    // linear operations on functions (the ranges and argument formatting must match)
    const asymptoticLF operator + (const asymptoticLF& rhs) const;
    const asymptoticLF operator - (const asymptoticLF& rhs) const;
    
	const asymptoticLF& operator += (const asymptoticLF& rhs);
	const asymptoticLF& operator -= (const asymptoticLF& rhs);

    // linear operations on all elements by a constant
    const asymptoticLF operator + (const valT& rhs) const;
    const asymptoticLF operator - (const valT& rhs) const;
    
    const asymptoticLF& operator += (const valT rhs);
    const asymptoticLF& operator -= (const valT rhs);

    // rescaling (element by element - including derivatives)
	const asymptoticLF operator * (const asymptoticLF& rhs) const;
	const asymptoticLF operator / (const asymptoticLF& rhs) const;
    const asymptoticLF& operator *= (const asymptoticLF& rhs);
    const asymptoticLF& operator /= (const asymptoticLF& rhs);
    

    // rescaling by a constant (including derivatives)
	const asymptoticLF operator * (const valT& rhs) const;
	const asymptoticLF operator / (const valT& rhs) const;
    const asymptoticLF& operator *= (const valT rhs);
    const asymptoticLF& operator /= (const valT rhs);
    
	// comparison (format and all elements)
	const bool operator == (const asymptoticLF& rhs);
	const bool operator != (const asymptoticLF& rhs);

    // reading and saving / file operations
    // binary formats for complete accuracy:
	int ReadF(ifstream*);
    int SaveF(ofstream*);
    
    // saving text format (for humans)
    // the saved file can be read by
    // ReadF(ifstream*);
    // however, the accuracy is lost
    int SaveFtxt(ofstream*);
    
    // reading and saving / file operations
    // binary formats for complete accuracy:
    int ReadF(string str)
    {
        ifstream fstr(str.c_str());
        return ReadF(&fstr);
    }
    // reading and saving / file operations
    // binary formats for complete accuracy:
    int ReadF(ifstream& istr)
    {
        return ReadF(&istr);
    }
 
    // reading and saving / file operations
    // binary formats for complete accuracy:
    int SaveF(string str)
    {
        ofstream fstr(str.c_str());
        return SaveF(&fstr);
    }
    // reading and saving / file operations
    // binary formats for complete accuracy:
    int SaveF(ofstream& istr)
    {
        return SaveF(&istr);
    }
    
    // saving text format (for humans)
    // the saved file can be read by
    // ReadF(ifstream*);
    // however, the accuracy is lost
    int SaveFtxt(string str)
    {
        ofstream fstr(str.c_str());
        return SaveFtxt(&fstr);
    }

    // saving text format (for humans)
    // the saved file can be read by
    // ReadF(ifstream*);
    // however, the accuracy is lost
    int SaveFtxt(ofstream& fstr)
    {
        return SaveFtxt(&fstr);
    }
    
    // direct return of the dataset
    // risky to ruin the allocation
    // should not modify/ delete, etc.
    storage<valT>* DirectAccess()
    {
        return interpolationF<valT>::values;
    }

};


// now two specific types
typedef asymptoticLF<float> asymptoticLF_float;
typedef asymptoticLF<double> asymptoticLF_double;


// this class defines unique asymptoticLF for complex type.
// instead of values being "complex"
// the values are used as real and this makes better for interpolations
// since the interpolation is not done in 2D, but is in 1D.
// otherwise it is more or less the same as asymptoticLF template
class asymptoticLF_complexv
{
  public:
    
    // construtors like of asymptoticLF
	asymptoticLF_complexv();
	asymptoticLF_complexv(const asymptoticLF_complexv& idat);
	asymptoticLF_complexv(storage<complexv>& idat,double ixmin,double istep);
	asymptoticLF_complexv(storage<complexv>& idat,double ixmin,double istep, complexv iderivI, complexv ivalF, complexv iderivF);

    // unique constructor
    // real values are stored in idatR
    // imaginary values are stored in idatI
    // ... and so on
    asymptoticLF_complexv(
                storage<double>& idatR,
                storage<double>& idatI,
                double ixmin,
                double istep,
                double iderivIr,
                double iderivIi,
                double ivalFr,
                double ivalFi,
                double iderivFr,
                double iderivFi
        );
    
    // destructor
	~asymptoticLF_complexv();


    // mathematical operations
    const asymptoticLF_complexv& operator = (const complexv& rhs);// this is meaningful
	const asymptoticLF_complexv& operator = (const asymptoticLF_complexv& rhs);
        
	const asymptoticLF_complexv& operator += (const asymptoticLF_complexv& rhs);
	const asymptoticLF_complexv& operator -= (const asymptoticLF_complexv& rhs);

	const asymptoticLF_complexv& operator += (const complexv& rhs);
	const asymptoticLF_complexv& operator -= (const complexv& rhs);
        
    // scaling
    const asymptoticLF_complexv& operator *= (const double& rhs);
	const asymptoticLF_complexv& operator /= (const double& rhs);

    const asymptoticLF_complexv& operator *= (const complexv& rhs);
	const asymptoticLF_complexv operator * (const complexv& rhs) const;

    // linear operations
	const asymptoticLF_complexv operator + (const asymptoticLF_complexv& rhs) const;
	const asymptoticLF_complexv operator - (const asymptoticLF_complexv& rhs) const;

        const asymptoticLF_complexv operator + (const complexv& rhs) const;
	const asymptoticLF_complexv operator - (const complexv& rhs) const;
	const asymptoticLF_complexv operator * (const double& rhs) const;
	const asymptoticLF_complexv operator / (const double& rhs) const;


	// comparison
	const bool operator == (const asymptoticLF_complexv& rhs);
	const bool operator != (const asymptoticLF_complexv& rhs);


    // file I/O operations are the same as asymptoticLF
	int ReadF(ifstream* fstr);
    int SaveF(ofstream* fstr);
    int SaveFtxt(ofstream* fstr);
    
    
    int SaveFtxt(ofstream& fstr)
    {
        return SaveFtxt(&fstr);
    }
    int SaveFtxt(string str)
    {
        ofstream fstr(str.c_str());
        return SaveFtxt(&fstr);
    }


    int SaveF(string str)
    {
        ofstream fstr(str.c_str());
        return SaveF(&fstr);
    }
    int SaveF(ofstream& fstr)
    {
        return SaveF(&fstr);
    }
    
    int ReadF(string str)
    {
        ifstream fstr(str.c_str());
        return ReadF(&fstr);
    }

	int ReadF(ifstream& fstr)
    {
        return ReadF(&fstr);
    }
    
    
    // returns the derivative of the function at "infinity"
    complexv GetD1Inf()
    {
        return complexv(RE->a,IM->a);
    }

    // sets causality options for real and imaginary parts separately
    // the causality is defined in the sense of asymptoticLF
    void SetCausality(int cr,int ci)
    {
        RE->causality = cr;
        IM->causality = ci;
    }

    // returns the number of explicitly given samples
    int GetN()
    {
        return RE->GetN();
    }
    
    // return the xmax value
    double GetXF()
    {
        return RE->GetXF();
    }

    // return the real part step value
    double GetStep()
    {
        return RE->GetStep();
    }

    
    //Returns value of the function for arbitrary argument
    complexv Get(double x)
    {
        return complexv(RE->Get(x),IM->Get(x));
    }
    
    complexv GetDiracAmplitude()
    {
             return complexv(RE->GetDiracAmplitude(),IM->GetDiracAmplitude());
    };
    complexv SetDiracAmplitude(complexv ival)
    {
        RE->SetDiracAmplitude(ival.real());
        RE->SetDiracAmplitude(ival.imag());
    };
    void SetDiracAmplitude(double ival)
    {
        RE->SetDiracAmplitude(ival);
    };


    // returns the asymptoticLF of real values
    asymptoticLF_double* DirectAccessRE()
    {
        return RE;
    }

    // returns the asymptoticLF of imaginary values
    asymptoticLF_double* DirectAccessIM()
    {
        return IM;
    }

    // returns the value of the function explicitly at xmin (no interpolation)
    complexv Get(string str)
    {
        return complexv(RE->Get(str),IM->Get(str));
    }
    
    
    void UpdateAxis(double imin, double istep)
    {
        RE->UpdateAxis(imin,istep);
        IM->UpdateAxis(imin,istep);
    }
    
  private:

    // the real and imaginary parts are stored in these functions
	asymptoticLF_double* RE;
	asymptoticLF_double* IM;
};

// Asymptotic linear function implementation
// for the above specified functions


//////////////////////////////
// now functions of templates

template<class valT> asymptoticLF<valT>::asymptoticLF()
:interpolationF<valT>()
{
    a=0;
    b=0;
    c=0;
    d=0;
    causality = 0;
    deltaamplitude = 0;
}
template<class valT> asymptoticLF<valT>::asymptoticLF(const asymptoticLF<valT>& rhs)
{
    int& num = interpolationF<valT>::num_sampl;
    
    num = rhs.values->GetSize();
    interpolationF<valT>::values = new storage<valT>(1);
    interpolationF<valT>::values->Allocate(num);
    
    interpolationF<valT>::values->SetBar(rhs.values->data1D,num);
    
    interpolationF<valT>::xmin = rhs.xmin;
    interpolationF<valT>::step = rhs.step;
    interpolationF<valT>::xmax = rhs.xmax;
    a=rhs.a;
    b=rhs.b;
    c=rhs.c;
    d=rhs.d;
    causality = rhs.causality;
    deltaamplitude = rhs.deltaamplitude;
    interpolationF<valT>::title = rhs.title;
    
}




template<class valT> asymptoticLF<valT>::asymptoticLF(valT* idat, double ixmin,double idx,int inump)
:interpolationF<valT>(idat,ixmin,idx,inump)
{
    int& num = interpolationF<valT>::num_sampl;
    double& step = interpolationF<valT>::step;
    
    //cout<<num;
    //num = idat.GetSize();
    if(num<2)
    {
        cout<<"Error: the data is not sufficient to create asymptoticLF\n";
    }
    else
    {
        //f=a*step+b
        
        a=(idat[num-1]-idat[num-2])/step;
        b=a*step+idat[num-1];
        d=idat[0];
        c=(idat[1]-idat[0])/step;
        
        causality = 0;
        deltaamplitude = 0;
    }
    
}


template<class valT> asymptoticLF<valT>::asymptoticLF(storage<valT>& idat,double ixmin,double idx)
:interpolationF<valT>(idat,ixmin,idx)
{
    int& num = interpolationF<valT>::num_sampl;
    double& step = interpolationF<valT>::step;
    
    //cout<<num;
    //num = idat.GetSize();
    if(num<2)
    {
        cout<<"Error: the data is not sufficient to create asymptoticLF\n";
    }
    else
    {
        //f=a*step+b
        
        a=(idat.data1D[num-1]-idat.data1D[num-2])/step;
        b=a*step+idat.data1D[num-1];
        d=idat.data1D[0];
        c=(idat.data1D[1]-idat.data1D[0])/step;
        
        causality = 0;
        deltaamplitude = 0;
        
    }
}
template<class valT> asymptoticLF<valT>::asymptoticLF(storage<valT>& idat,double ixmin,double idx, valT iderivI, valT ivalF, valT iderivF)
:interpolationF<valT>(idat,ixmin,idx)
{
    a=iderivF;
    b=ivalF;
    d=idat.data1D[0];
    c=iderivI;
    causality = 0;
        deltaamplitude = 0;
    
}

template<class valT> asymptoticLF<valT>::~asymptoticLF()
{
    a=0;
    b=0;
    c=0;
    d=0;
    causality = 0;
        deltaamplitude = 0;
}


template<class valT>
valT asymptoticLF<valT>::Get(double x)
{
    // there is a bug in this function:
    // when argument is equal to xmax, the function returns nonsence.
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //cout<<"THERE IS AN ERROR!!!!!\n";
    // fixed (?) more testing is needed
    //
    // update:
    // the bug seems to be solved
    
    
    if(causality == 3 && x<0)
    {
        // odd functionality for real functions f(-x) = -f(x)
        return -Get(-x);
    }
    
    if(causality == 2 && x<0)
    {
        // even functionality for real functions f(-x) = f(x)
        return Get(-x);
    }
    
    
    
    int& num = interpolationF<valT>::num_sampl;
    double& step = interpolationF<valT>::step;
    double& xmin = interpolationF<valT>::xmin;
    double& xmax = interpolationF<valT>::xmax;
    
    
    
    if(step<0)
    {
        double dist = xmin - x;
        
        asymptoticLF<valT> tfunct(*interpolationF<valT>::values,xmin,-step, -c, b, -a);
        return tfunct.Get(xmin+dist);
    }
    
    //    if(x==  14000.0+((16000.0-14000.0)/1000)*999 )
    //    {
    //        cout<<x<<"\n";
    //
    //  cout<<xmin<<"\n";
    //  cout<<xmax<<"\n";
    //  cout<<step<<"\n";
    //  cout<<num<<"\n";
    //}
    double eps = step / 1.0e10;
    
    if( x < xmin+eps )
        return c*(x-xmin)+d;
    
    //std::cout.precision(15);
    //std::cout<<"x from asymptoticLF: "<<fixed<<x<<"\n";
    
    if(x > xmax-eps ) //!!!!!!! this fix may be temporary
        //if(x>=(xmax-step) ) //!!!!!!! this fix may be temporary
        return a*(x-xmax)+b;
    
    
    
    // getting index
    double interval = xmax-xmin;
    if(interval<eps)
        cout<<"Strange interval\n";
    double z=(x-xmin)/interval*num;
    
    int ii= (int)floor(z);
    
    // getting displacement from ii
    z=z-(double)ii;
    
    //cout<<"ii: "<<ii<<"\n";
    //cout<<"z: "<<z<<"\n";
    
    
    if(fabs(z) <eps )
        // no interpolation is necessary
        return interpolationF<valT>::values->data1D[ii];
    
    valT v0;// = values->data1D[i-1];
    valT v1;// = values->data1D[i];
    valT v2;// = values->data1D[i+1];
    valT v3;// = values->data1D[i+2];
    // Bezier type interpolation between points v1 and v2
    // 1 finding all necessary derivatives:
    valT deriva;// = v1-v0;
    valT derivb;// = v2-v1;
    valT derivc;// = v3-v2;
    valT deriv1;// = (deriva+derivb)/2.0;
    valT deriv2;// = (derivc+derivb)/2.0;
    
    
    if(num == 1)
    {
        
        v1 = interpolationF<valT>::values->data1D[0];
        v2 = b;
        deriv1 = a*step;
        deriv2 = c*step;
        
    }
    
    
    
    else if(num == 2)
    {
        
        v1 = interpolationF<valT>::values->data1D[0];
        v2 = interpolationF<valT>::values->data1D[1];
        deriv1 = a*step;
        deriv2 = c*step;
        
        
    }
    else
    {
        // here three or more
        
        
        
        if(ii==0)
        {
            //v0 = values->data1D[0];
            v1 = interpolationF<valT>::values->data1D[0];
            v2 = interpolationF<valT>::values->data1D[1];
            v3 = interpolationF<valT>::values->data1D[2];
            derivb = v2-v1;
            deriva = c*step; //derivb;
            derivc = v3-v2;
            
            deriv1 = deriva;//(deriva+derivb)/2.0; //
            deriv2 = (derivc+derivb)/2.0;
        }
        else if(ii==num-1)
        {
            v0 = interpolationF<valT>::values->data1D[num-2];
            v1 = interpolationF<valT>::values->data1D[num-1];
            v2 = b;
            //v3 = values->data1D[ii+1];
            derivb = v2-v1;
            deriva = v1-v0;
            derivc = a*step; //derivb
            
            deriv1 = (deriva+derivb)/2.0;
            deriv2 = derivc;//(derivc+derivb)/2.0;
        }
        
        
        else if(ii==num-2)
        {
            v0 = interpolationF<valT>::values->data1D[num-3];
            v1 = interpolationF<valT>::values->data1D[num-2];
            v2 = interpolationF<valT>::values->data1D[num-1];
            v3 = b;
            derivb = v2-v1;
            deriva = v1-v0;
            derivc = v3-v2;
            
            deriv1 = (deriva+derivb)/2.0;
            deriv2 = (derivc+derivb)/2.0;
            
        }
        else if(ii>0 && ii<num-2)
        {
            // when four points are available
            v0 = interpolationF<valT>::values->data1D[ii-1];
            v1 = interpolationF<valT>::values->data1D[ii];
            v2 = interpolationF<valT>::values->data1D[ii+1];
            v3 = interpolationF<valT>::values->data1D[ii+2];
            derivb = v2-v1;
            deriva = v1-v0;
            derivc = v3-v2;
            
            deriv1 = (deriva+derivb)/2.0;
            deriv2 = (derivc+derivb)/2.0;
            
            //cout<<"here\n";
            
        }
    }
    
    // Bezier type interpolation between points v1 and v2
    // 1 finding all necessary derivatives:
    // finding values p1 and p2
    valT p1 = v1+deriv1*0.333333333333333;
    valT p2 = v2-deriv2*0.333333333333333;
    
    //if(ii==num-1)
    //        {
    //
    //
    //        cout<<"\n"<<v1<<"\n";
    //        cout<<"\n"<<v2<<"\n";
    //        cout<<"\n"<<p1<<"\n";
    //        cout<<"\n"<<p2<<"\n";
    //        }
    
    return  interpolationF<valT>::retBez(z, v1, p1, p2, v2);
    
}

template<class valT>
int asymptoticLF<valT>::SaveFtxt(ofstream* fstr)
{
    *fstr<<"# asymptoticLF<valT>, v. 3.0:\n";
    fstr->precision(15);
    
    interpolationF<valT>::SaveFtxt(fstr);
    
    *fstr<<"# Derivative at xmax:\n";
    *fstr<<a<<"\n";
    *fstr<<"# Value at xmax:\n";
    *fstr<<b<<"\n";
    *fstr<<"# Derivative at xmin:\n";
    *fstr<<c<<"\n";
    *fstr<<"# Value at xmin:\n";
    *fstr<<d<<"\n";
    *fstr<<"# Symmetry parameter ***0 -> ignored,  2 -> f(-t) = f(t), 3 -> f(-t) = -f(t) ***:\n";
    *fstr<<causality<<"\n";
    *fstr<<"# Amplitude of the dirac delta at zero:\n";
    *fstr<<deltaamplitude<<"\n";
    *fstr<<"# End of asymptoticLF<valT>, v. 3.0\n";
    
    return 0;
}

template<class valT>
int asymptoticLF<valT>::SaveF(ofstream* fstr)
{
    char message[100];
    
    memset ( message, 0, 100*sizeof(char) );
    char str1[]="% asymptoticLF<valT>, v. 3.0: bin\n";
    strcpy (message,str1);
    fstr->write(message, 100*sizeof(char));
    
    interpolationF<valT>::SaveF(fstr);
    fstr->write( (char*) &a,sizeof(valT));
    fstr->write( (char*) &b,sizeof(valT));
    fstr->write( (char*) &c,sizeof(valT));
    fstr->write( (char*) &d,sizeof(valT));
    fstr->write( (char*) &causality,sizeof(int));
    fstr->write( (char*) &deltaamplitude,sizeof(valT));
    
    
    memset ( message, 0, 100*sizeof(char) );
    char str2[]="\n% End of asymptoticLF<valT>, v. 3.0: bin\n";
    strcpy (message,str2);
    fstr->write(message, 100*sizeof(char));
    
    return 0;
}
template<class valT>
int asymptoticLF<valT>::ReadF(ifstream* fstr)
{
    if (!fstr->is_open())
    {
		cout<<"Error: asymptoticLF<valT>: File has not been opened. Aborting.\n";
		return 2;
    }

        

	bool binary;
	char ccc = fstr->get();

	  if ( (ccc >= '0') && (ccc <= '9') )
	  {
	    fstr->putback (ccc);
	    // a number
	    binary = false;
	  }
	  else if (ccc == '#' )
	  {
		  // a comment
		  fstr->putback (ccc);
		  binary = false;
	  }
      else if (ccc == '\n' )
      {
          // a comment
          fstr->putback (ccc);
          binary = false;
      }
      else if (ccc == '\t' )
      {
          // a comment
          fstr->putback (ccc);
          binary = false;
      }
	  else if (ccc == '%' )
	  	  {
	  		  // a binary
	  		  fstr->putback (ccc);
	  		  binary = true;
	  	  }
	  	  else
	  {
			// something else
			fstr->putback (ccc);
		    binary = false;
		    cout<<"Warning: interpolationF<valT>::ReadF : unrecognized initial character; trying to read as text\n";
	  }


	    if(binary)
	    {
	    	int numsize;
	        char cmessage[100];
	        //fstr->getline(cmessage, 100);
	        // reading full message
	        fstr->read( cmessage,100*sizeof(char));
	        string message(cmessage);

	        if(message==string("% asymptoticLF<valT>, v. 3.0: bin\n"))
	        {
	        	cout<<"# Reading "<<message;
        
//        // reading message
//        char cmessage[100];
//        fstr->getline(cmessage, 100);
//        string message(cmessage);
//        cout<<"Reading "<<message<<"\n";
//
//        if(message == string("# Asymptotic Linear Function, v. 2.0: bin"))
//        {
//
//            // reading the remaining zeros
//            char str1[]="# Asymptotic Linear Function, v. 2.0: bin\n";
//            int len1 = sizeof(str1);
//            //cout<<len1;
//            fstr->read( cmessage,(100-len1+1)*sizeof(char));
            
            
	        	// inner interpolation function
	        	interpolationF<valT>::ReadF(fstr);
            
	        	//cout<<"I am here\n";
            
	        	fstr->read( (char*) &a,sizeof(valT));
	        	fstr->read( (char*) &b,sizeof(valT));
	        	fstr->read( (char*) &c,sizeof(valT));
	        	fstr->read( (char*) &d,sizeof(valT));
	        	fstr->read( (char*) &causality,sizeof(int));
	        	fstr->read( (char*) &deltaamplitude,sizeof(valT));
	        	char Nmessage[100];
	        	fstr->read(Nmessage, 100*sizeof(char));

	        	//cout<<"I am here\n";

	        	// successfull reading
            	return 0;
	        }
            else
            {
        		cout<<"Error: asymptoticLF<valT>: Wrong file format. Aborting.\n";
        		return 2;
            }

        }

	    else
	    {
	    	// reading text format (comments can be included)
	        // reading full text file



        //else if (message == string("# Asymptotic Linear Function, v. 2.0: text:short"))
        //{
        //
        //    // inner interpolation function
        //    interpolationF<valT>::ReadF(fstr);
        //    a=0;
        //    b=0;
        //    c=0;
        //    d=0;
        //    causality = 0;
        //}
        //else if (message == string("# Asymptotic Linear Function, v. 2.0: text:full"))
        //{
            
            toolsIO tio;
            string str = "";
            
            // inner interpolation function
            interpolationF<valT>::ReadF(fstr);
            
            tio.StreamSkipTrailers(fstr);
            getline(*fstr,str);
            tio.StreamSkipTrailers(str);
            a = tio.fromString<double>(str);
            
            tio.StreamSkipTrailers(fstr);
            getline(*fstr,str);
            tio.StreamSkipTrailers(str);
            b = tio.fromString<double>(str);
            
            tio.StreamSkipTrailers(fstr);
            getline(*fstr,str);
            tio.StreamSkipTrailers(str);
            c = tio.fromString<double>(str);
            
            tio.StreamSkipTrailers(fstr);
            getline(*fstr,str);
            tio.StreamSkipTrailers(str);
            d = tio.fromString<double>(str);
            
            
            tio.StreamSkipTrailers(fstr);
            getline(*fstr,str);
            tio.StreamSkipTrailers(str);
            causality = tio.fromString<int>(str);
            
            tio.StreamSkipTrailers(fstr);
            getline(*fstr,str);
            tio.StreamSkipTrailers(str);
            deltaamplitude = tio.fromString<double>(str);
            // successfull reading
            return 0;

        }
        //else
        //{
        //    cout<<"Error: asymptoticLF<valT>::ReadF(ifstream* fstr): file format incorrect.\n";
        //}
        
        
        
    //}
    //else
    //{
    //    cout<<"Error: asymptoticLF<valT>::ReadF(ifstream* fstr): file not open for reading\n";
    //}
    
    return 0;
}


//---------------

template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator = (const asymptoticLF<valT>& rhs)
{
    if (this != &rhs)
    {
        
        
        interpolationF<valT>::xmax=rhs.xmax;
        interpolationF<valT>::xmin=rhs.xmin;
        interpolationF<valT>::step=rhs.step;
        
        if(interpolationF<valT>::num_sampl != rhs.num_sampl)
        {
            interpolationF<valT>::values->Delete();
            interpolationF<valT>::values->Allocate(rhs.num_sampl);
            interpolationF<valT>::num_sampl = rhs.num_sampl;
        }
        interpolationF<valT>::values->SetBar(rhs.values->data1D,interpolationF<valT>::num_sampl);
        
        a = rhs.a;
        b = rhs.b;
        c = rhs.c;
        d = rhs.d;
        causality = rhs.causality;
        deltaamplitude = rhs.deltaamplitude;
    }
    return *this;    // Return ref for multiple assignment
}



template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator = ( interpolationF<valT>& rhs)
{
    if (true)
    {


        interpolationF<valT>::xmax=rhs.GetXF();
        interpolationF<valT>::xmin=rhs.GetXI();
        interpolationF<valT>::step=rhs.GetStep();

        if(interpolationF<valT>::num_sampl != rhs.GetN())
        {
            interpolationF<valT>::num_sampl = rhs.GetN();
            interpolationF<valT>::values->Delete();
            interpolationF<valT>::values->Allocate(interpolationF<valT>::num_sampl);
        }
        valT* dataset = rhs.DirectAccessD();
        interpolationF<valT>::values->SetBar(dataset,interpolationF<valT>::num_sampl);

        a = 0.0;
        b = 0.0;
        c = 0.0;
        d = 0.0;
        causality = 0.0;
        deltaamplitude = 0.0;
    }
    return *this;    // Return ref for multiple assignment
}



template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator = (const valT rhs)
{
    interpolationF<valT>::operator=(rhs);
    a = 0.0;
    b = rhs;
    c = 0.0;
    d = rhs;
    causality = 2;
    deltaamplitude = 0;
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator += (const valT rhs)
{
    interpolationF<valT>::operator+=(rhs);
    //a += rhs;
    b += rhs;
    //c += rhs;
    d += rhs;
    
    return *this;
}
template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator -= (const valT rhs)
{
    interpolationF<valT>::operator-=(rhs);
    //a -= rhs;
    b -= rhs;
    //c -= rhs;
    d -= rhs;
    
    return *this;
}
template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator *= (const valT rhs)
{
    interpolationF<valT>::operator*=(rhs);
    a *= rhs;
    b *= rhs;
    c *= rhs;
    d *= rhs;
    
    deltaamplitude *= rhs;
    
    return *this;
}
template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator /= (const valT rhs)
{
    interpolationF<valT>::operator/=(rhs);
    a /= rhs;
    b /= rhs;
    c /= rhs;
    d /= rhs;

    deltaamplitude /= rhs;
    
    return *this;
}



template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator += (const asymptoticLF<valT>& rhs)
{
    interpolationF<valT>::operator+=(rhs);
    a += rhs.a;
    b += rhs.b;
    c += rhs.c;
    d += rhs.d;
    
    return *this;
}
template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator -= (const asymptoticLF<valT>& rhs)
{
    interpolationF<valT>::operator-=(rhs);
    a -= rhs.a;
    b -= rhs.b;
    c -= rhs.c;
    d -= rhs.d;
    
    return *this;
    
}
template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator *= (const asymptoticLF<valT>& rhs)
{
    cout<<"Error: operator* cannot be used in asymptoticLF<valT>\n";
    
    return *this;
}
template<class valT>
const asymptoticLF<valT>& asymptoticLF<valT>::operator /= (const asymptoticLF<valT>& rhs)
{
    cout<<"Error: operator/ cannot be used in asymptoticLF<valT>\n";
    
    return *this;
    
}





template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator * (const asymptoticLF<valT>& rhs)
const
{
    cout<<"Error: operator* cannot be used in asymptoticLF<valT>\n";
    
    return asymptoticLF<valT>(*this);
}

template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator / (const asymptoticLF<valT>& rhs)
const
{
    cout<<"Error: operator/ cannot be used in asymptoticLF<valT>\n";
    
    return asymptoticLF<valT>(*this);
}

template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator + (const asymptoticLF<valT>& rhs)
const
{
    return asymptoticLF<valT>(*this) += rhs;
}

template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator - (const asymptoticLF<valT>& rhs)
const
{
    return asymptoticLF<valT>(*this) -= rhs;
}

template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator + (const valT& rhs)
const
{
    return asymptoticLF<valT>(*this) += rhs;
}

template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator - (const valT& rhs)
const
{
    return asymptoticLF<valT>(*this) -= rhs;
}

template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator * (const valT& rhs)
const
{
    return asymptoticLF<valT>(*this) *= rhs;
}

template<class valT>
const asymptoticLF<valT> asymptoticLF<valT>::operator / (const valT& rhs)
const
{
    return asymptoticLF<valT>(*this) /= rhs;
}


template<class valT>
const bool asymptoticLF<valT>::operator == (const asymptoticLF<valT>& rhs)
{
    if(this == &rhs)
        return true;
    
    bool result = true;
    
    result = result && interpolationF<valT>::num_sampl == rhs.num_sampl;
    result = result && interpolationF<valT>::xmax == rhs.xmax;
    result = result && interpolationF<valT>::xmin == rhs.xmin;
    result = result && interpolationF<valT>::step == rhs.step;
    
    if(result == false) return false;
    
    for(int ind=0; ind<interpolationF<valT>::num_sampl; ind ++)
        result = result && (interpolationF<valT>::values->data1D[ind]==rhs.values->data1D[ind]);
    
    if(result == false) return false;
    
    result = result && (a == rhs.a);
    result = result && (b == rhs.b);
    result = result && (c == rhs.c);
    result = result && (d == rhs.d);
    result = result && (causality == rhs.causality);
    result = result && (deltaamplitude == rhs.deltaamplitude);
    
    return result;    // Return ref for multiple assignment
}

template<class valT>
const bool asymptoticLF<valT>::operator != (const asymptoticLF<valT>& rhs)
{
    return !(*this==rhs);    // Return ref for multiple assignment
}

