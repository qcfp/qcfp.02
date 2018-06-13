#include<iostream>
#include <time.h>

#include"calculator_3rd_secular_cumulant.hpp"

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

    
#ifdef MPIPROCESSING
    
    int ierr, num_procs, my_id;
    ierr = MPI_Init(&argc, &argv);
    
    /* find out MY process ID, and how many processes were started. */
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    if(my_id == 0)
        printf("Hello! %i processes has been started in MPI\n", num_procs);
    
#endif
    
    
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

    // perform calculations based on filename:
    calculator_3rd_secular_cumulant calculator(ifilename);

    // finally run calculation
    calculator.LaunchBasic();

#ifdef MPIPROCESSING
    if(my_id == 0)
    {
#endif
  
    // print the result
    calculator.Publish(ofilename);

    
//    // additional printout
//    if(true)
//    {
//        string name;
//        ofstream ofs;
//        
//        name = ofilename + ".2DRe.txt";
//        ofs.open(name.c_str());
//        calculator.Publish2DRe(ofs);
//        ofs.close();
//
//
//        name = ofilename + ".2DIm.txt";
//        ofs.open(name.c_str());
//        calculator.Publish2DIm(ofs);
//        ofs.close();
//    }
        
#ifdef MPIPROCESSING
}
#endif

        
//    time(&time2);
//    cout << "Execution time was " << difftime(time2, time1) << " seconds.\n";
    
#ifdef MPIPROCESSING
    // printout will be done only from a single master process
    ierr = MPI_Finalize();
#endif
    return 0;
}
