#pragma once

#include<cmath>
#include <typeinfo>       // operator typeid


#include<iostream>
#include<fstream>
#include<string>

using namespace std;



#include"../storage/storage.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../complexv/complexv.hpp"


// Interpolated function based on the given dataset:
// main part is stored in the array $values of arbitrary size
// Argument in this part goes from 0 to $xmax and the data is sampled using num_sampl points.
// After $xmax the function (argument is $x)  is zero.

template <class valT> 
class interpolationF
{
    friend class asymptoticLF_complexv;
    
    
  public:

    // constructors:

    // creates an empty function which cannot provide values
	interpolationF();

	// creates a function from the given
	interpolationF(const interpolationF&);

	// creates a function from the given dataset with provided xmin and step
	interpolationF(storage<valT>& idat,double ix,double dx);

	// creates a function from the given array with provided xmin and step
	interpolationF(valT* idat,double ix,double dx,int nump);

	// creates a zero function with provided xmin, step and number of points
	interpolationF(double ix,double dx,int nump);

	// destroys the function
	~interpolationF();

	//Returns value of the function for arbitrary argument
	valT Get(double x);

	//operators
	// assignment
	const interpolationF& operator = (valT rhs);// this is meaningful 
	//interpolationF& operator = (int rhs);// this is meaningful only when rhs = 0
	//interpolationF& operator = (double rhs);// this is meaningful only when rhs = 0
	const interpolationF& operator = (const interpolationF& rhs);
        
	const interpolationF& operator += (const valT rhs);
	const interpolationF& operator -= (const valT rhs);
	const interpolationF& operator *= (const valT rhs);
	const interpolationF& operator /= (const valT rhs);
	const interpolationF& operator += (const interpolationF& rhs);
	const interpolationF& operator -= (const interpolationF& rhs);
	const interpolationF& operator *= (const interpolationF& rhs);
	const interpolationF& operator /= (const interpolationF& rhs);


	// multiplication by a scalar
	const interpolationF operator * (const valT& rhs) const;
	const interpolationF operator / (const valT& rhs) const;
	const interpolationF operator + (const valT& rhs) const;
	const interpolationF operator - (const valT& rhs) const;

	// linear operations
	const interpolationF operator + (const interpolationF& rhs) const;
	const interpolationF operator - (const interpolationF& rhs) const;
	const interpolationF operator * (const interpolationF& rhs) const;
	const interpolationF operator / (const interpolationF& rhs) const;


	// comparison
	const bool operator == (const interpolationF& rhs);
	const bool operator != (const interpolationF& rhs);
    
		// saving function to a file in a binary format
		int SaveF(string);

		// reading function from a file in a binary or text format
		int ReadF(string);

		// saving function to a file in a binary format
		int SaveF(ofstream*);

		// reading function from a file in a binary or text format
		int ReadF(ifstream*);

		// saving function to a file in a binary format
		int SaveF(ofstream&);

		// reading function from a file in a binary or text format
		int ReadF(ifstream&);
    
		// saving function to a file in a binary format
		int SaveFtxt(ofstream*);

		// saving function to a file in a binary format
    	int SaveFtxt(ofstream&);

    	// saving function to a file in a binary format
    	int SaveFtxt(string);

        // sets the stored values to be equal to the ones stored in irep
    	// be careful not to mess up the memory
    	// the size of given values must be not larger than the storage space
    	void DirectAssign(interpolationF<valT>& irep);

        // replaces the stored array of data
        // be careful not to mess up the memory
        valT* DirectExchange(valT* irep);
        
        // returns the stored array of data
        // be careful not to mess up the memory
        valT* DirectAccessD();
        
        // changes xmin and xmax and step values
        void UpdateAxis(double imin, double istep)
        {
            xmin = imin;
            step = istep;
            xmax = xmin +num_sampl*step;
        }

        // Returns the number of points
        int GetN()
        {
            return num_sampl;
        }

        // returns the initial argument xmin
        double GetXI()
        {
            return xmin;
        }

        // returns the final argument xmax
        double GetXF()
        {
            return xmax;
        }

        // returns the argument stepsize
        double GetStep()
        {
        	return step;
        }

        void SetTitle(std::string istr)
        {
            title = istr;
        }
        std::string GetTitle()
        {
            return title;
        }

        bool CompareType(const interpolationF& rhs);



  protected:
      	storage<valT>* values;

        int num_sampl;

        //step size of argument during initial part
        double xmax;
        double xmin;
        double step;
        
        std::string title;
        
  //private:
        // bezier interpolation
        valT  retBez(double&, valT&, valT&, valT&, valT&);


};



#pragma once

//////////////////////////////
// now functions of templates

template<class valT>
interpolationF<valT>::interpolationF()
{
    num_sampl = 0;
    values = new storage<valT>(1);
    
    xmax = 0.0;
    xmin = 0.0;
    step = 0.0;
    title = "";
}

template<class valT>
interpolationF<valT>::interpolationF(storage<valT>& idat,double ix,double dx)
{
    num_sampl = idat.GetSize();
    values = new storage<valT>(1);
    values->Allocate(num_sampl);
    values->SetBar(idat.data1D,num_sampl);
    
    xmin = ix;
    step = dx;
    xmax = ix+step*num_sampl;
    title = "";
}
template<class valT>
interpolationF<valT>::interpolationF(double ix,double dx,int nump)
{
    num_sampl = nump;
    values = new storage<valT>(1);
    values->Allocate(num_sampl);
    
    xmin = ix;
    step = dx;
    xmax = ix+step*nump;
    title = "";
}
template<class valT>
interpolationF<valT>::interpolationF(valT* idat,double ix,double dx,int nump)
{
    num_sampl = nump;
    values = new storage<valT>(1);
    values->Allocate(num_sampl);
    values->SetBar(idat,num_sampl);
    
    xmin = ix;
    step = dx;
    xmax = ix+step*nump;
    title = "";
}
template<class valT>
interpolationF<valT>::interpolationF(const interpolationF<valT>& rhs)
{
    num_sampl = rhs.values->GetSize();
    values = new storage<valT>(1);
    values->Allocate(num_sampl);
    
    values->SetBar(rhs.values->data1D,num_sampl);
    
    xmin = rhs.xmin;
    step = rhs.step;
    xmax = rhs.xmax;
    title = rhs.title;

}
template<class valT>
interpolationF<valT>::~interpolationF()
{
    values->Delete();
    delete values;
    num_sampl = 0;
    xmax = 0.0;
    xmin = 0.0;
    step = 0.0;
    title = "";
}

template<class valT>
bool interpolationF<valT>::CompareType(const interpolationF<valT>& rhs)
{
    
    bool condition= this != &rhs;
    if(condition)
    {
        condition = condition && this->xmax==rhs.xmax;
        condition = condition &&  this->xmin==rhs.xmin;
        condition = condition &&  this->step==rhs.step;
        condition = condition &&  this->num_sampl == rhs.num_sampl;
    }
    else return true;
    
    return condition;
}


template<class valT>
valT interpolationF<valT>::Get(double x)
{
    
    //    cout<<xmin<<"\n";
    //    cout<<step<<"\n";
    //    cout<<xmax<<"\n";
    
    if(step<0)
    {
        double dist = xmin - x;
        
        interpolationF<valT> tfunct(*values,xmin,-step);
        return tfunct.Get(xmin+dist);
    }
    
    
    if(num_sampl == 1)
    {
        return values->data1D[0];
    }
    
    double eps = step / 1.0e10;
    
    if(x<xmin+eps)
        return values->data1D[0];
    
    if(x>xmax- step-eps)
        return values->data1D[num_sampl-1];
    
    
    
    // getting index
    double interval = xmax-xmin;
    double z=(x-xmin)/interval*num_sampl;
    
    int ii= (int)floor(z);
    
    // getting displacement from ii
    z=z-(double)ii;
    
    
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
    
    
    if(num_sampl == 2)
    {
        // linear interpolation
        //    valT a1,b1;
        //     b1=values->data1D[0];
        //	 a1=values->data1D[1]-values->data1D[0];
        //	return a1*z+b1;
        
        v1 = values->data1D[0];
        v2 = values->data1D[1];
        deriv1 = 0.0;
        deriv2 = 0.0;
        
        
    }
    else
    {
        // here three or more
        if(ii>0 && ii<num_sampl-2)
        {
            // when four points are available
            v0 = values->data1D[ii-1];
            v1 = values->data1D[ii];
            v2 = values->data1D[ii+1];
            v3 = values->data1D[ii+2];
            deriva = v1-v0;
            derivb = v2-v1;
            derivc = v3-v2;
            deriv1 = (deriva+derivb)/2.0;
            deriv2 = (derivc+derivb)/2.0;
            
        }
        else if(ii==0)
        {
            v0 = values->data1D[0];
            v1 = values->data1D[0];
            v2 = values->data1D[1];
            v3 = values->data1D[2];
            derivb = v2-v1;
            derivc = v3-v2;
            deriv1 = 0.0;
            deriv2 = (derivc+derivb)/2.0;
        }
        else if(ii==num_sampl-2)
        {
            v0 = values->data1D[ii-1];
            v1 = values->data1D[ii];
            v2 = values->data1D[ii+1];
            v3 = values->data1D[ii+1];
            deriva = v1-v0;
            derivb = v2-v1;
            deriv1 = (deriva+derivb)/2.0;
            deriv2 = 0.0;
            
        }
        
    }
    
    // Bezier type interpolation between points v1 and v2
    // 1 finding all necessary derivatives:
    // finding values p1 and p2
    valT p1 = v1+deriv1*0.333333333333333;
    valT p2 = v2-deriv2*0.333333333333333;
    
    x=z;
    
    return  retBez(x, v1, p1, p2, v2);
    
}

template<class valT>
valT  interpolationF<valT>::retBez(double& x, valT& a1, valT& a2, valT& a3, valT& a4)
{
    return  (1-x)*(1-x)*(1-x)*a1+3*(1-x)*(1-x)*x*a2+3*(1-x)*x*x*a3+x*x*x*a4;
}

template<class valT>
int interpolationF<valT>::SaveF(ofstream& fstr)
{
    return SaveF(&fstr);
}
template<class valT>
int interpolationF<valT>::SaveF(ofstream* fstr)
{
    char message[100];
    memset ( message, 0, 100*sizeof(char) );
    char str1[]="% interpolationF<valT>, v. 3.0: bin\n";
    strcpy (message,str1);
    fstr->write(message, 100*sizeof(char));


    const char* cstr = title.c_str();
    int cstrl = title.size();
    //cout<<"title = "<<title<<" "<<"size: "<<cstrl<<"\n";
    fstr->write( (char*) &cstrl,sizeof(int));
    fstr->write( cstr,cstrl*sizeof(char));



    fstr->write( (char*) &num_sampl,sizeof(int));
    fstr->write( (char*) &xmin,sizeof(double));
    //fstr->write( (char*) &xmax,sizeof(double));
    fstr->write( (char*) &step,   sizeof(double));
    
    if(num_sampl>0)
        fstr->write( (char*) values->data1D,num_sampl*sizeof(valT));
    
    //xmax = step*num_sampl;
    xmax = num_sampl*step+xmin;

    memset ( message, 0, 100*sizeof(char) );
    char str2[]="\n% End of interpolationF<valT>, v. 3.0: bin\n";
    strcpy (message,str2);
    fstr->write(message, 100*sizeof(char));
    return 0;
}
template<class valT>
int interpolationF<valT>::SaveFtxt(ofstream* fstr)
{
    fstr->precision(15);
    
    *fstr<<"# interpolationF<valT>, v. 3.0:\n";
    *fstr<<"# Number of samples:\n";
    *fstr<<num_sampl<<"\n";
    *fstr<<"# Initial argument xmin:\n";
    *fstr<<xmin<<"\n";
    //*fstr<<"# Final argument xmax:\n";
    //*fstr<<xmax<<"\n";
    *fstr<<"# Argument step size:\n";
    *fstr<<step<<"\n";
    
    *fstr<<"# The values of the function at the predefined points (N values):\n";

    // for double, float, int... other are undefined
    	// error is not reported
        if(num_sampl>0)
        {
            for(int ind = 0; ind<num_sampl; ind++)
                *fstr<< values->data1D[ind]<<"\n";
        }
    *fstr<<"# End of interpolationF<valT>, v. 3.0:\n";
    return 0;
}



template<class valT>
int interpolationF<valT>::SaveFtxt(ofstream& fstr)
{
    return SaveFtxt(&fstr);
}
template<class valT>
int interpolationF<valT>::SaveFtxt(string str)
{
    ofstream fstr(str.c_str());
    return SaveFtxt(&fstr);
}


template<class valT>
int interpolationF<valT>::SaveF(string str)
{
    
    ofstream fstr(str.c_str());
    return SaveF(&fstr);
}
template<class valT>
int interpolationF<valT>::ReadF(ifstream& fstr)
{
    return ReadF(&fstr);
}
template<class valT>
int interpolationF<valT>::ReadF(string str)
{
    ifstream fstr(str.c_str());
    return ReadF(fstr);
}
template<class valT>
int interpolationF<valT>::ReadF(ifstream* fstr)
{
    
    if (!fstr->is_open())
    {
		cout<<"Error: interpolationF<valT>: File has not been opened. Aborting.\n";
		return 2;
    }


	bool binary;
	char ccc = fstr->get();

	  if ( (ccc >= '0') && (ccc <= '9') )
	  {
	    fstr->putback (ccc);
	    // a number: will be read
	    binary = false;
	  }
	  else if (ccc == '#' )
	  {
		  // a comment: will be ignored
		  fstr->putback (ccc);
		  binary = false;
	  }
      else if (ccc == '\n' )
      {
          // end of line character: will be removed
          fstr->putback (ccc);
          binary = false;
      }
      else if (ccc == '\t' )
      {
          // end of line character: will be removed
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
		    cout<<"Warning: interpolationF<valT>::ReadF : unexpected initial character; trying to read as text\n";
	  }



    
    if(binary)
    {
    	int numsize;
        char cmessage[100];
        //fstr->getline(cmessage, 100);
        // reading full message
        fstr->read( cmessage,100*sizeof(char));
        string message(cmessage);

        if(message==string("% interpolationF<valT>, v. 3.0: bin\n"))
        {
        	cout<<"# Reading "<<message;
        	// reading the remaining zeros
        	//char str1[]="# Interpolation Function, v. 2.0: bin\n";
        	//int len1 = sizeof(str1);
        	//fstr->read( cmessage,(100-len1+1)*sizeof(char));


                int cstrl;
            	fstr->read( (char*) &cstrl,sizeof(int));
                char* cstr = new char[cstrl+1];
                memset ( cstr, 0, (cstrl+1)*sizeof(char) );
                fstr->read( cstr,cstrl*sizeof(char));

                const char* ccstr = cstr;
                title = string(ccstr);
                
    //cout<<"title = "<<title<<"\n";
                
        	fstr->read( (char*) &numsize,sizeof(int));
        	fstr->read( (char*) &xmin,sizeof(double));
        		//fstr->read( (char*) &xmax,sizeof(double));
        	fstr->read( (char*) &step,sizeof(double));
        	num_sampl = numsize;
            xmax = num_sampl*step+xmin;
    		values->Delete();
        	if(numsize>0)
        	{
        		//cout<<numsize<<"\n";
        		values->Allocate(num_sampl);
        		fstr->read( (char*) values->data1D,num_sampl*sizeof(valT));
        	}
        	else
        	{
        		cout<<"Warning: interpolationF<valT>: Trying to read an empty function. Aborting.\n";
        		//return 1;
        	}
        
        	// here reading the ending message
        	char Nmessage[100];
        	fstr->read(Nmessage, 100*sizeof(char));
                
                delete[] cstr;

        	// successfull reading
        	return 0;
        }
        else
        {
    		cout<<"Error: interpolationF<valT>: Wrong file format. Aborting.\n";
    		return 2;
        }
    }
    else
    {
    	// reading text format (comments can be included)
        // reading short text file
        
        toolsIO tio;
        string str = "";
        
        // reading the number of values in the function
        tio.StreamSkipTrailers(fstr);
        getline(*fstr,str);
        tio.StreamSkipTrailers(str);
        num_sampl = tio.fromString<int>(str);
        
        // reading the initial argument value xmin
        tio.StreamSkipTrailers(fstr);
        getline(*fstr,str);
        tio.StreamSkipTrailers(str);
        xmin = tio.fromString<double>(str);

        // reading the physical step between the values
        tio.StreamSkipTrailers(fstr);
        getline(*fstr,str);
        tio.StreamSkipTrailers(str);
        step = tio.fromString<double>(str);
        
        // double
        	// reading all values from the file
        	//storage<double> arr(2); // two dimensional
        	//arr.Allocate(num_sampl,1);
        	//if(tio.ReadRectangular(&arr, *fstr))
        	//{
        	//	cout<<"Error: interpolationF<valT>: Some problem with text file reading\n";
        	//	return 2;// 1;
        	//}
        	//// saving
        	//values->Delete();
        	//values->Allocate(num_sampl);
        	//for(int iii=0; iii<num_sampl; iii++)
        	//	values->data1D[iii]=arr.data2D[iii][0];

        values->Delete();
        if(num_sampl>0)
        {
            values->Allocate(num_sampl);
            for(int ind = 0; ind<num_sampl; ind++)
            {
                tio.StreamSkipTrailers(fstr);
                *fstr>> values->data1D[ind];
            }
        }

        	//xmin = 0;
        xmax = num_sampl*step +xmin;

        return 0;

    }
    //else
    //    cout<<"Error: Cannot read this format\n";
    //
    return 2;
}


template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator = (const interpolationF<valT>& rhs)
{
    if (this != &rhs)
    {
        xmax=rhs.xmax;
        xmin=rhs.xmin;
        step=rhs.step;
        
        if(num_sampl != rhs.num_sampl)
        {
            //cout<<num_sampl<<" numsample\n";
            values->Delete();
            num_sampl = rhs.num_sampl;
            values->Allocate(num_sampl);
            //cout<<num_sampl<<" numsample\n";
        }
        values->SetBar(rhs.values->data1D,num_sampl);
        
    }
    return *this;    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator = (const valT rhs)
{
    {
        xmax=1;
        xmin=0;
        step=1;
        
        values->Delete();
        values->Allocate(1);
        num_sampl = 1;
        
        // getting lhs data object
        values->data1D[0] = rhs;
        
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator += (const valT rhs)
{
    {
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]+=rhs;
        
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator -= (const valT rhs)
{
    {
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]-=rhs;
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator *= (const valT rhs)
{
    {
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]*=rhs;
        
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator /= (const valT rhs)
{
    {
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]/=rhs;
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator += (const interpolationF<valT>& rhs)
{
    if(CompareType(rhs))
    {
        
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]+=rhs.values->data1D[ind];
    }
    else
    {
        cout<<"Error: cannot append object (file interpolationF.hpp)\n";
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator -= (const interpolationF<valT>& rhs)
{
    if(CompareType(rhs))
    {
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]-=rhs.values->data1D[ind];
    }
    else
    {
        cout<<"Error: cannot append object (file interpolationF.hpp)\n";
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator *= (const interpolationF<valT>& rhs)
{
    if(CompareType(rhs))
    {
        
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]*=rhs.values->data1D[ind];
    }
    else
    {
        cout<<"Error: cannot append object (file interpolationF.hpp)\n";
    }
    return *this;    // Return ref for multiple assignment
}
template<class valT>
const interpolationF<valT>& interpolationF<valT>::operator /= (const interpolationF<valT>& rhs)
{
    if(CompareType(rhs))
    {
        for(int ind=0; ind<num_sampl; ind++)
            this->values->data1D[ind]/=rhs.values->data1D[ind];
    }
    else
    {
        cout<<"Error: cannot append object (file interpolationF.hpp)\n";
    }
    return *this;    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT> interpolationF<valT>::operator + (const interpolationF<valT>& rhs)
const
{
    return interpolationF<valT>(*this += rhs);    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT> interpolationF<valT>::operator - (const interpolationF<valT>& rhs)
const
{
    return interpolationF<valT>(*this -= rhs);    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT> interpolationF<valT>::operator * (const interpolationF<valT>& rhs)
const
{
    return interpolationF<valT>(*this *= rhs);    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT> interpolationF<valT>::operator / (const interpolationF<valT>& rhs)
const
{
    return interpolationF<valT>(*this /= rhs);    // Return ref for multiple assignment
}


template<class valT>
const interpolationF<valT> interpolationF<valT>::operator + (const valT& rhs)
const
{
    return interpolationF<valT>(*this += rhs);    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT> interpolationF<valT>::operator - (const valT& rhs)
const
{
    return interpolationF<valT>(*this -= rhs);    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT> interpolationF<valT>::operator * (const valT& rhs)
const
{
    return interpolationF<valT>(*this *= rhs);    // Return ref for multiple assignment
}

template<class valT>
const interpolationF<valT> interpolationF<valT>::operator / (const valT& rhs)
const
{
    return interpolationF<valT>(*this /= rhs);    // Return ref for multiple assignment
}


template<class valT>
const bool interpolationF<valT>::operator == (const interpolationF<valT>& rhs)
{
    if(this == &rhs)
        return true;
    
    bool result = CompareType(rhs);
    
    if(result)
    {
        
        for(int ind=0; ind<this->num_sampl; ind ++)
            result = result && (this->values->data1D[ind]==rhs.values->data1D[ind]);
        
    }
    return result;    // Return ref for multiple assignment
}

template<class valT>
const bool interpolationF<valT>::operator != (const interpolationF<valT>& rhs)
{
    return !(operator ==(rhs));    // Return ref for multiple assignment
}


template<class valT>
void interpolationF<valT>::DirectAssign(interpolationF<valT>& irep)
{
    if(num_sampl<irep.num_sampl)
    {
        cout<<"Error: in void interpolationF<valT>::DirectReplace(interpolationF<valT>& irep)\n";
        return;
    }
    else
    {
        values->SetBar(irep.values->data1D,irep.num_sampl);
    }
}
template<class valT>
valT* interpolationF<valT>::DirectExchange(valT* irep)
{
    valT* tmp = values->data1D;
    values->data1D = irep;
    return tmp;
}


template<class valT>
valT* interpolationF<valT>::DirectAccessD()
{
    return values->data1D;
}

