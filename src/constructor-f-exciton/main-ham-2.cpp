// this program reads a triangular one-exciton hamiltonian block
// from a file
// constructs the double-exciton hamiltonian block
// assuming two-level sites
// calculates eigenstates and eigenvectors of both matrices
// and writes then data into  files

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
    
    // for two-level molecules
    int numel2 = numel*(numel-1)/2;

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
    
    storage<double>* ham2 = new storage<double>(2); // two dimensional
	ham2->Allocate(numel2,numel2);
    obj.GetHam2(ham2);
//    for(int ii=0;   ii<numel2; ii++)
//        for(int ij=0;   ij<numel2; ij++)
//            cout<<ii<<" "<<ij<<" "<<ham2->data2D[ii][ij]<<"\n";
    

    string fham2 = filename + ".ham2";
    tio.WriteTriangular(ham2,fham2);
    ham2->Delete();
    delete ham2;

    storage<double>* evals2 = new storage<double>(1); // one dimensional
    evals2->Allocate(numel2);
    obj.GetEvals2(evals2);
    string fevals2 = filename + ".evals2";
    tio.WriteRectangular(evals2,fevals2);
    evals2->Delete();
    delete evals2;
    
    storage<double>* evecs2 = new storage<double>(2); // two dimensional
	evecs2->Allocate(numel2,numel2);
    obj.GetEvecs2(evecs2);
    string fevecs2 = filename + ".evecs2";
    tio.WriteRectangular(evecs2,fevecs2);
    evecs2->Delete();
    delete evecs2;

    
    // that is it.
}

