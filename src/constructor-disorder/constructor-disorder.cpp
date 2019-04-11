// everything is performed
#include"constructor-disorder.hpp"

constructor_disorder::constructor_disorder(int iset)
{
    sval = iset;
    processor = new toolsRandom(sval);
}
constructor_disorder::constructor_disorder()
{
    sval = 0;
    processor = new toolsRandom();
}

constructor_disorder::~constructor_disorder()
{
    delete processor;
}

int constructor_disorder::DiagIntra(storage<double>& ham, double sigma)
{
    int numl,numr;
    
    ham.GetSize(numl,numr);

    if(numl!=numr || numl==0 || numr == 0 )
    {
        cout<<"Error: wrong size of the matrix\n";
        return 1;
    }
    
    for(int iv = 0; iv<numl; iv++)
    {
        ham.data2D[iv][iv] += processor->RandGD(0.0, sigma);
    }
    return 0;
}


int constructor_disorder::DiagInter(storage<double>& ham, double sigma)
{
    int numl,numr;
    
    ham.GetSize(numl,numr);

    if(numl!=numr || numl==0 || numr == 0 )
    {
        cout<<"Error: wrong size of the matrix\n";
        return 1;
    }
    
    double dval = processor->RandGD(0.0, sigma);
    
    for(int iv = 0; iv<numl; iv++)
    {
        ham.data2D[iv][iv] += dval;
    }
    return 0;
}


