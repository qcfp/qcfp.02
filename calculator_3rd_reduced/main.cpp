/* 
 * File:   main.cpp
 * Author: dariusa
 *
 * Created on November 6, 2015, 8:18 PM
 */

#include <cstdlib>

#include"../propagatorExciton/propagatorExciton.hpp"
#include"calculator_3rd_reduced.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    
    
    if (!argc ==3 )
    {
        cout<<"Error: calling sequence is \n";
        cout<<"program ifile ofile\n";
        return 0;
    }
    
    
    toolsIO tio;

    string ifilename(argv[1]);
    string ofilename(argv[2]);

//    string ifilename("input.txt");
//    string ofilename("sig.txt");

        calculator_3rd_reduced<propagatorExciton> calculator(ifilename);

        // finally run calculation
        calculator.LaunchBasic();

        // print the result
        calculator.Publish(ofilename);

    return 0;
}

