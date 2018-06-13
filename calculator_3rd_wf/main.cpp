#include<iostream>
#include <time.h>

#include"../propagatorExcitonWF/propagatorExcitonWF.hpp"
#include"calculator_3rd_wf.hpp"

using std::cin;
using std::cout;



#include"../toolsIO/toolsIO.hpp"

#include<string>
using std::string;

#ifdef MPIPROCESSING
#include<mpi.h>
#endif



int main(int argc,  char * argv[])
{


    
    
    //    time_t time1, time2;
//    time(&time1);
    if(argc != 3 )
    {    cout<<"Error: specify the input file and the output file\n";
        return 0;
    }

    toolsIO tio;

    string ifilename(argv[1]);
    string ofilename(argv[2]);

//    string ifilename("input.txt");
//    string ofilename("sig.txt");
    ifstream istr(ifilename.c_str());
    if(!istr.is_open())
    {
        cout<<"Error: input file cannot be read\nquitting.\n";
        return 0;
    }



    // perform calculations based on filename:
    calculator_3rd_wf<propagatorExcitonWF,wavefunctionD1> calculator(ifilename);


    // finally run calculation
    calculator.LaunchBasic();

    // print the result
    calculator.Publish(ofilename);

    
//    // additional printout
//    if(true)
//    {
//        string name;
//        ofstream ofs;
        
//        name = ofilename + ".2DRe.txt";
//        ofs.open(name.c_str());
//        calculator.Publish2DRe(ofs);
 //       ofs.close();


 //       name = ofilename + ".2DIm.txt";
 //       ofs.open(name.c_str());
 //       calculator.Publish2DIm(ofs);
 //       ofs.close();
 //   }

//    time(&time2);
//    cout << "Execution time was " << difftime(time2, time1) << " seconds.\n";
    

    return 0;
}
