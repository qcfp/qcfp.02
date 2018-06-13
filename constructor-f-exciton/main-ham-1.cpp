// this program reads a triangular one-exciton hamiltonian block
// from a file
// calculates eigenstates and eigenvectors 
// and writes then data into five files

#include"constructor-f-exciton.hpp"

#include"../eigen/eigen.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../storage/storage.hpp"

#include<string>
using std::string;


int main(int argc, const char * argv[])
{

    // reads the matrix and the number of elements in the matrix
    if(argc != 3 )
    {    
        cout<<"Error: specify the filename for a triangular hamiltonian matrix and the number of sites\n";
        return 1;
    }
 
    toolsIO tio;
    
    string filename(argv[1]);
    int numel = tio.fromString<int>(string(argv[2]));
    
    if(numel<=0)
    {
        cout<<"Error: The rank of the matrix must be positive integer\n";
        return 1;
    }
    
    storage<double>* arr = new storage<double>(2); // two dimensional
	arr->Allocate(numel,numel);

    if(tio.ReadTriangular(arr, filename))
    {
        cout<<"Error: Some problem with file reading\n";
        return 1;
    }
    
    // adding the other triangle:
    for(int indl=0; indl<numel; indl++)
    {
        for(int indr=0; indr<indl; indr++)
        {
            arr->data2D[indr][indl] = arr->data2D[indl][indr] ;
        }
    }
    
       
       
    ////////
    constructor_f_exciton obj(arr);
    arr->Delete();
    delete arr;
    ////////
    
    storage<double>* evals = new storage<double>(1); // one dimensional
    evals->Allocate(numel);
    obj.GetEvals(evals);
    string fevals = filename + ".evals";
    tio.WriteRectangular(evals,fevals);
    evals->Delete();
    delete evals;

    storage<double>* evecs = new storage<double>(2); // two dimensional
	evecs->Allocate(numel,numel);
    obj.GetEvecs(evecs);
    string fevecs = filename + ".evecs";
    tio.WriteRectangular(evecs,fevecs);
    evecs->Delete();
    delete evecs;
    
    // that is it.
}

