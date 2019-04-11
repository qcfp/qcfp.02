#include"complexv.hpp"
#include<iostream>
#include<string>

using std::cout;
using std::string;

void complexv::printb(ostream* str)
{
    
    double re = this->real();
    double im = this->imag();
    
//    str->write (label, 9);
    str->write ((char*)&re, sizeof(re));
    str->write ((char*)&im, sizeof(im));
}
void complexv::printbn(ostream* str)
{
    
    double re = this->real();
    double im = this->imag();
    
    //str->write (label, 9);
    str->write ((char*)&re, sizeof(re));
    str->write ((char*)&im, sizeof(im));
}
void complexv::readbn(istream* str)
{
    
    double re;// = this->real();
    double im;// = this->imag();
    
    str->read ((char*)&re, sizeof(re));
    str->read ((char*)&im, sizeof(im));
    (*this)=complexv(re,im);
}
void complexv::readb(istream* str)
{
    double re;// = this->real();
    double im;// = this->imag();
    //string strtest("qcfptcomx");
    //char istr[10];
    
    //str->read (istr,9);
    
    //if(string(istr)==strtest)
    {   
        str->read ((char*)&re, sizeof(re));
        str->read ((char*)&im, sizeof(im));
        (*this)=complexv(re,im);
    }
    //else
    //{
    //    cout<<"Error: cannot read complex data\n";
    //}
}


  complexv complexv::conjugate( complexv& z)
{
    return complexv(z.real(),-z.imag());
}
complexv complexv::conjugate()
{
    return complexv(this->real(),-this->imag());
}

