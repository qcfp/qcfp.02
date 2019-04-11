// functions for complex variable


#include"asymptoticLF.hpp"


asymptoticLF_complexv::asymptoticLF_complexv()
{
    RE = new asymptoticLF_double;
    IM = new asymptoticLF_double;
}
asymptoticLF_complexv::asymptoticLF_complexv(const asymptoticLF_complexv& idat)
{
    int num = idat.RE->num_sampl;
    if(num>0)
    {
        RE = new asymptoticLF_double;
        IM = new asymptoticLF_double;
        *RE = *(idat.RE);
        *IM = *(idat.IM);
    }
    else
    {
        RE = new asymptoticLF_double;
        IM = new asymptoticLF_double;
    }
}


asymptoticLF_complexv::asymptoticLF_complexv(storage<complexv>& idat,double ixmin,double idx)
{
    int num = idat.GetSize();
    
    if(num>0)
    {
        
        storage<double> vR(1);
        storage<double> vI(1);
        vR.Allocate(num);
        vI.Allocate(num);
        
        for(int ind = 0; ind<num; ind++)
        {
            vR.data1D[ind] = idat.data1D[ind].real();
            vI.data1D[ind] = idat.data1D[ind].imag();
        }
        
        RE = new asymptoticLF_double(vR,ixmin,idx);
        IM = new asymptoticLF_double(vI,ixmin,idx);
    }
    else
    {
        RE = new asymptoticLF_double;
        IM = new asymptoticLF_double;
    }
    
}

asymptoticLF_complexv::asymptoticLF_complexv(
                                             storage<complexv>& idat,
                                             double ixmin,
                                             double idx,
                                             complexv iderivI,
                                             complexv ivalF,
                                             complexv iderivF)
{
    int num = idat.GetSize();
    
    if(num>0)
    {
        storage<double> vR(1);
        storage<double> vI(1);
        vR.Allocate(num);
        vI.Allocate(num);
        
        for(int ind = 0; ind<num; ind++)
        {
            vR.data1D[ind] = idat.data1D[ind].real();
            vI.data1D[ind] = idat.data1D[ind].imag();
        }
        
        RE = new asymptoticLF_double(vR,ixmin,idx);
        IM = new asymptoticLF_double(vI,ixmin,idx);
        
        RE->a = iderivF.real();
        RE->b = ivalF.real();
        RE->c = iderivI.real();
        RE->d = idat.data1D[0].real();
        
        IM->a = iderivF.imag();
        IM->b = ivalF.imag();
        IM->c = iderivI.imag();
        IM->d = idat.data1D[0].imag();
        
    }
}
asymptoticLF_complexv::asymptoticLF_complexv(
                                             storage<double>& idatR,
                                             storage<double>& idatI,
                                             double ixmin,
                                             double idx,
                                             double iderivIr,
                                             double iderivIi,
                                             double ivalFr,
                                             double ivalFi,
                                             double iderivFr,
                                             double iderivFi
                                             )
{
    int num = idatR.GetSize();
    
    if(num>0)
    {
        RE = new asymptoticLF_double(idatR,ixmin,idx,iderivIr, ivalFr, iderivFr);
        IM = new asymptoticLF_double(idatI,ixmin,idx,iderivIi, ivalFi, iderivFi);
        
    }
    else
    {
        RE = new asymptoticLF_double;
        IM = new asymptoticLF_double;
    }
    
}
asymptoticLF_complexv::~asymptoticLF_complexv()
{
    //
    //cout<<"test inside asymptoticLF_complexv\n";
    delete RE;
    delete IM;
}


//------------------

const asymptoticLF_complexv& asymptoticLF_complexv::operator = (const complexv& rhs)
{
    *RE = rhs.real();
    *IM = rhs.imag();
    return *this;
}
const asymptoticLF_complexv& asymptoticLF_complexv::operator = (const asymptoticLF_complexv& rhs)
{
    if(*this != rhs)
    {
        *RE = *(rhs.RE);
        *IM = *(rhs.IM);
    }
    return *this;
}




const asymptoticLF_complexv& asymptoticLF_complexv::operator += (const asymptoticLF_complexv& rhs)
{
    *RE += *(rhs.RE);
    *IM += *(rhs.IM);
    return *this;
}

const asymptoticLF_complexv& asymptoticLF_complexv::operator -= (const asymptoticLF_complexv& rhs)
{
    *RE -= *(rhs.RE);
    *IM -= *(rhs.IM);
    return *this;
}


const asymptoticLF_complexv& asymptoticLF_complexv::operator += (const complexv& rhs)
{
    *RE += rhs.real();
    *IM += rhs.imag();
    return *this;
}

const asymptoticLF_complexv& asymptoticLF_complexv::operator -= (const complexv& rhs)
{
    *RE -= rhs.real();
    *IM -= rhs.imag();
    return *this;
}
const asymptoticLF_complexv& asymptoticLF_complexv::operator *= (const double& rhs)
{
    *RE *= rhs;
    *IM *= rhs;
    return *this;
}

const asymptoticLF_complexv& asymptoticLF_complexv::operator *= (const complexv& rhs)
{
    asymptoticLF_double rtmp = *RE;
    *RE = rtmp*rhs.real()-(*IM)*rhs.imag();
    *IM = rtmp*rhs.imag()+(*IM)*rhs.real();
    return *this;
}

const asymptoticLF_complexv& asymptoticLF_complexv::operator /= (const double& rhs)
{
    *RE /= rhs;
    *IM /= rhs;
    return *this;
}







const asymptoticLF_complexv asymptoticLF_complexv::operator + (const asymptoticLF_complexv& rhs)
const
{
    return asymptoticLF_complexv(*this) += rhs;
}

const asymptoticLF_complexv asymptoticLF_complexv::operator - (const asymptoticLF_complexv& rhs)
const
{
    return asymptoticLF_complexv(*this) -= rhs;
}


const asymptoticLF_complexv asymptoticLF_complexv::operator + (const complexv& rhs)
const
{
    return asymptoticLF_complexv(*this) += rhs;
}

const asymptoticLF_complexv asymptoticLF_complexv::operator - (const complexv& rhs)
const
{
    return asymptoticLF_complexv(*this) -= rhs;
}
const asymptoticLF_complexv asymptoticLF_complexv::operator * (const double& rhs)
const
{
    return asymptoticLF_complexv(*this) *= rhs;
}

const asymptoticLF_complexv asymptoticLF_complexv::operator * (const complexv& rhs)
const
{
    return asymptoticLF_complexv(*this) *= rhs;
}

const asymptoticLF_complexv asymptoticLF_complexv::operator / (const double& rhs)
const
{
    return asymptoticLF_complexv(*this) /= rhs;
}





// comparison
const bool asymptoticLF_complexv::operator == (const asymptoticLF_complexv& rhs)
{
    return (*RE == *(rhs.RE)) && (*IM == *(rhs.IM));
}
const bool asymptoticLF_complexv::operator != (const asymptoticLF_complexv& rhs)
{
    return ! (*this == rhs);
}

int asymptoticLF_complexv::SaveFtxt(ofstream* fstr)
{
    *fstr<<"# asymptoticLF_complexv, v. 3.0:\n";
    fstr->precision(15);
    
    RE->SaveFtxt(fstr);
    IM->SaveFtxt(fstr);
    
    *fstr<<"# End of asymptoticLF_complexv, v. 3.0:\n";
    
    return 0;
}


int asymptoticLF_complexv::SaveF(ofstream* fstr)
{
    char message[100];// = "# Complex asymptotic function, v. 2.0\n";
    
    memset ( message, 0, 100*sizeof(char) );
    char str1[]="% asymptoticLF_complexv, v. 3.0: bin\n";
    strcpy (message,str1);
    fstr->write(message, 100*sizeof(char));
    
    RE->SaveF(fstr);
    IM->SaveF(fstr);
    
    memset ( message, 0, 100*sizeof(char) );
    char str2[]="% End of asymptoticLF_complexv, v. 3.0: bin\n";
    strcpy (message,str2);
    fstr->write(message, 100*sizeof(char));
    
    return 0;
}
int asymptoticLF_complexv::ReadF(ifstream* fstr)
{
    if (!fstr->is_open())
    {
		cout<<"Error: asymptoticLF_complexv: File has not been opened. Aborting.\n";
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
      else if (ccc == '\t' )
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
    		char message[100];// = "# Complex asymptotic function, v. 2.0\n";
    		fstr->read(message, 100*sizeof(char));
    		cout<<"# Reading "<<string(message);
    		RE->ReadF(fstr);
    		IM->ReadF(fstr);
    		char Nmessage[100];// = "# End of Complex asymptotic function, v. 2.0\n";
    		fstr->read(Nmessage, 100*sizeof(char));
    	}
    	else
    	{
    		cout<<"# Reading asymptoticLF_complexv text format\n";
    		RE->ReadF(fstr);
    		IM->ReadF(fstr);
    	}


    return 0;
}
// void asymptoticLF_complexv::SetCausality(int ic)
// {
//     RE->causality = ic;
// }
