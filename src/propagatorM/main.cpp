#include<iostream>

#include"propagatorM.hpp"
#include"../toolsIO/toolsIO.hpp"

using std::cin;
using std::cout;
using std::ifstream;
using std::ofstream;


// this programs solves the master equation

// first the program reads the number of values (dimensionality of the problem)
// then it reads a matrix as transfer rates
// then program reads a vector as initial values,
// then it reads the output timestep
// and finally it reads the number of points on the output

// then it writes the new values after propagation into the output file


int main(int argc, const char * argv[])
{
    toolsIO tio;

    if(argc != 3 )
    {    cout<<"Error: specify the input file and the output file\n";
        return 0;
    }

    string ifilename(argv[1]);
    string ofilename(argv[2]);

    storage<double> matr(2);

    

    ifstream ifs(ifilename.c_str());

    string str = "";

    cout<<"Reading the rank of the transfer tensor\n";
    // read the size of the matrix
    tio.StreamSkipTrailers(&ifs);
    getline(ifs,str);
    tio.StreamSkipTrailers(str);
    int numEl = tio.fromString<int>(str);

    matr.Allocate(numEl,numEl);

    // The very transfer tensor 
    cout<<"Reading the transfer tensor\n";
    tio.StreamSkipTrailers(&ifs);
    if( tio.ReadRectangular(&matr, ifs) )
    {
        cout<<"Error: Some problem with file reading\n";
        return 0;// 1;
    }

    // reading initial condition
    cout<<"Reading the initial vector\n";
    storage<double> vec1(1);
    vec1.Allocate(numEl);

    tio.StreamSkipTrailers(&ifs);
    if( tio.ReadRectangular(&vec1, ifs) )
    {
        cout<<"Error: Some problem with file reading\n";
        return 0;// 1;
    }


    // reading the output time step:
    cout<<"Reading the output timestep\n";
    storage<double> stimestep(1); 
    stimestep.Allocate(1);
    tio.StreamSkipTrailers(&ifs);
    if( tio.ReadRectangular(stimestep, ifs) )
    {
        cout<<"Error: Some problem with file reading\n";
        return 0;// 1;
    }

    cout<<"Reading the output number of points\n";
    // reading the output number of iterations:
    storage<int> snumber(1); 
    snumber.Allocate(1);
    tio.StreamSkipTrailers(&ifs);
    if( tio.ReadRectangular(snumber, ifs) )
    {
        cout<<"Error: Some problem with file reading\n";
        return 0;// 1;
    }

    
    ifs.close();

    
    std::cout<<"Calculating\n";


    propagatorM obj;
    
    //obj.InitMaster(num, matr);
    //obj.PropagateMaster(vec1,numt,dt);
    //obj.ExitMaster();
    
    obj.InitMaster(matr);
    ofstream ofs(ofilename.c_str());
    
    for(int it=0;it<snumber.data1D[0];it++){
        double time = stimestep.data1D[0]*it;
        storage<double> result = obj.Get(vec1,time);
        ofs<<"#--------------------- at time = "<<time<<"\n";
        // printing result to output file:
        for(int ind1 = 0; ind1<numEl; ind1 ++ )
            ofs<<result.data1D[ind1]<<"\n";
    }   

    
    
}
