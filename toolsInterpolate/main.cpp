#include<iostream>

#include"toolsInterpolate.hpp"
#include<string>

using std::cout;


// this program reads a vector from standard input and prints 
// the interpolated value at the output
// four values are printed:
// 1 rational interpolation
// 2 error of the rational interpolation
// 3 polynomial interpolation
// 4 error of the polynomial interpolation


// first value is the number of points in the function
// second, third, ... are x and y values of the dataset
// last is the requested x value to be calculated

int main()
{

    std::cout<<"Enter number of points\n";

    int num;
    std::cin>>num;
    double* vec1; // for x values
    double* vec2; // for y values
    
    if(num>0)
    {
        vec1=new double[num];
        vec2=new double[num];
    }
    else
    {
        std::cout<<"Error: the number of points must be positive\n";
        return 0;
    }

    std::cout<<"Reading "<<num<<" points as x y pairs\n";
    
    // reading input
	for(int ind = 0; ind<num; ind ++ )
        std::cin>>vec1[ind] >> vec2[ind];
    
    
    // the requested point
    double xval;
    std::cin>>xval;
    double retv;
    double reterr;
    

    // doing interpolation
    toolsInterpolate lobj;
    
    lobj.ratint(vec1,vec2,num,xval,retv,reterr);
	std::cout<<retv<<"\t"<<reterr<<"\n";

    lobj.polint(vec1,vec2,num,xval,retv,reterr);
	std::cout<<retv<<"\t"<<reterr<<"\n";

    ////////
    // making test dataset
//    for(int ind=0; ind<1000; ind++)
//    {
//        double x=ind*1.0/100;
//        lobj.barint(vec1,vec2,num,x,retv,reterr);
//        std::cout<<x<<"\t"<<retv<<"\t"<<reterr<<"\n";
//    }
    

    // cleaning
    delete[] vec1;
    delete[] vec2;
    

}
