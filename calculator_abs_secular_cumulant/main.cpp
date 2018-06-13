#include<iostream>

#include"calculator_abs_secular_cumulant.hpp"

using std::cin;
using std::cout;

int main(int argc, const char * argv[])
{
    
    if(argc != 3 )
    {    cout<<"Error: specify the input file and the output file\n";
        return 0;
    }
    
    toolsIO tio;
    
    string ifilename(argv[1]);
    string ofilename(argv[2]);
    
    ifstream istr(ifilename.c_str());
    if(!istr.is_open())
    {
        cout<<"Error: input file cannot be read\nquitting.\n";
        return 0;
    }

    //    string ifilename("input.txt");
    //    string ofilename("sig.txt");
    
    // perform calculations based on filename:
    calculator_abs_secular_cumulant calculator(ifilename);
    
    // finally run calculation
    calculator.Launch();
    
    // print the result
    calculator.Publish(ofilename);
    
    return 0;
}
