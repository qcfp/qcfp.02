// this program reads a triangular matrix from a file
// calculates its eigenstates and eigenvectors
// and writes then into two files

#include"eigen.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../storage/storage.hpp"

#include<string>
using std::string;


int main(int argc, const char * argv[])
{

    if(argc != 3 )
    {    
        cout<<"Error: specify the filename for a triangular matrix and the rank of the matrix\n";
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
    
    storage<double>* evals = new storage<double>(1); // one dimensional
    evals->Allocate(numel);
       
    storage<double>* evecs = new storage<double>(2); // two dimensional
	evecs->Allocate(numel,numel);
       
    eigen obj;
    obj.GetEigenSolution(arr, evals, evecs);
    
    
    arr->Delete();
    delete arr;
 
    
    // write data to file
    string fevals = filename + ".evals";
    string fevecs = filename + ".evecs";
    
    tio.WriteRectangular(evals,fevals);
    evals->Delete();
    delete evals;

    tio.WriteRectangular(evecs,fevecs);
    evecs->Delete();
    delete evecs;

}

