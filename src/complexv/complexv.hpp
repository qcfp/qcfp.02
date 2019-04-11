#pragma once

// this class defines a complex number.
// it is characterized by a label 
// qcfptcomx
// qcfp--type--4labelcode

#include<complex>
using std::complex;
#include<iostream>
using std::ostream;
using std::istream;

//typedef complex<double> complexdouble;

class complexv: public complex<double>
{
// here some new improvements could be done

public:
    
    // constructors
    complexv() : complex<double>()     {
//        label[0]='q';
//        label[1]='c';
//        label[2]='f';
//        label[3]='p';
//        label[4]='t';
//        label[5]='c';
//        label[6]='o';
//        label[7]='m';
//        label[8]='x';
    }
    complexv(int x, int y) : complex<double>(x,y)     {
//        label[0]='q';
//        label[1]='c';
//        label[2]='f';
//        label[3]='p';
//        label[4]='t';
//        label[5]='c';
//        label[6]='o';
//        label[7]='m';
//        label[8]='x';
    }
    complexv(double x, double y) : complex<double>(x,y)    {
//        label[0]='q';
//        label[1]='c';
//        label[2]='f';
//        label[3]='p';
//        label[4]='t';
//        label[5]='c';
//        label[6]='o';
//        label[7]='m';
//        label[8]='x';
    }

    complexv(double x) : complex<double>(x)    {
//        label[0]='q';
//        label[1]='c';
//        label[2]='f';
//        label[3]='p';
//        label[4]='t';
//        label[5]='c';
//        label[6]='o';
//        label[7]='m';
//        label[8]='x';
    }

    complexv(complex<double> x) : complex<double>(x)    {
//        label[0]='q';
//        label[1]='c';
//        label[2]='f';
//        label[3]='p';
//        label[4]='t';
//        label[5]='c';
//        label[6]='o';
//        label[7]='m';
//        label[8]='x';
    }

    complexv conjugate( complexv&);
    complexv conjugate();
    double re(){
	return this->real();
	};
    double im(){
	return this->imag();
	};
    double norm(){
	double re = this->real();
	double im = this->imag();
	return re*re+im*im;
	};
    double abs(){
	return sqrt(norm());
    };
    // prints binary
    void printb(ostream*);
    // prints binary without label
    void printbn(ostream*);
    
    // reads binary
    void readb(istream*);
    // reads binary without label
    void readbn(istream*);


private:
    //char label[10];//={'q','c','f','p','t','c','o','m','x'};

    
};

// constant value
static const complexv coni(0,1);
static const complexv cnni(0,-1);


