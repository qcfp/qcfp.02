// this program reads a triangular one-exciton hamiltonian block
// from a file
// constructs the double-exciton hamiltonian block
// assuming two-level sites
// calculates eigenstates and eigenvectors of both matrices
// and writes then data into  files

#include"constructor-f-exciton.hpp"

#include"../eigen/eigen.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../storage/storage.hpp"

#include<string>
using std::string;


int main(int argc, const char * argv[])
{

    // reads the hamiltonian and dipole matrices and the number of sites
    if(argc != 4 )
    {    
        cout<<"Error: specify the filename for a triangular hamiltonian matrix, filename for a three-column matrix for the transition dipoles and the number o sites\n";
        return 1;
    }
 
    toolsIO tio;
    
    string filenameh(argv[1]);
    string filenamed(argv[2]);
    int numel = tio.fromString<int>(string(argv[3]));
    
    if(numel<=0)
    {
        cout<<"Error: The rank of the matrix must be positive integer\n";
        return 1;
    }
    
    // for two-level molecules
    int numel2 = numel*(numel-1)/2;

    storage<double>* arr = new storage<double>(2); // two dimensional
	arr->Allocate(numel,numel);

    storage<double>* ard = new storage<double>(2); // two dimensional
	ard->Allocate(numel,3);

    if(tio.ReadTriangular(arr, filenameh))
    {
        cout<<"Error: Some problem with file reading\n";
        return 1;
    }
    
    // adding the other triangle:
    for(int indl=0; indl<numel; indl++)
        for(int indr=0; indr<indl; indr++)
            arr->data2D[indr][indl] = arr->data2D[indl][indr];
    
    if(tio.ReadRectangular(ard, filenamed))
    {
        cout<<"Error: Some problem with file reading\n";
        return 1;
    }
    
    // making dipoles
    storage<dvector3d>* dip = new storage<dvector3d>(1); // two dimensional
	dip->Allocate(numel);
    for(int ind=0; ind<numel; ind ++)
        dip->data1D[ind]=dvector3d(ard->data2D[ind][0],ard->data2D[ind][1],ard->data2D[ind][2]);
    ard->Delete();
    delete ard;

    ////////
    constructor_f_exciton obj(arr,dip);
    arr->Delete();
    delete arr;
    dip->Delete();
    delete dip;
    ////////
    
    string name;
    
    storage<double>* evals = new storage<double>(1); // one dimensional
    evals->Allocate(numel);
    obj.GetEvals(evals);
    name = filenameh + ".evals";
    tio.WriteRectangular(evals,name);
    evals->Delete();
    delete evals;

    storage<double>* evecs = new storage<double>(2); // two dimensional
	evecs->Allocate(numel,numel);
    obj.GetEvecs(evecs);
    name = filenameh + ".evecs";
    tio.WriteRectangular(evecs,name);
    evecs->Delete();
    delete evecs;
    
    storage<double>* ham2 = new storage<double>(2); // two dimensional
	ham2->Allocate(numel2,numel2);
    obj.GetHam2(ham2);
    name = filenameh + ".ham2";
    tio.WriteTriangular(ham2,name);
    ham2->Delete();
    delete ham2;

    storage<double>* evals2 = new storage<double>(1); // one dimensional
    evals2->Allocate(numel2);
    obj.GetEvals2(evals2);
    name = filenameh + ".evals2";
    tio.WriteRectangular(evals2,name);
    evals2->Delete();
    delete evals2;
    
    storage<double>* evecs2 = new storage<double>(2); // two dimensional
	evecs2->Allocate(numel2,numel2);
    obj.GetEvecs2(evecs2);
    name = filenameh + ".evecs2";
    tio.WriteRectangular(evecs2,name);
    evecs2->Delete();
    delete evecs2;

    
    // transition dipoles single
    storage<dvector3d>* di = new storage<dvector3d>(1); // one dimensional
	di->Allocate(numel);
    storage<double>* ar1 = new storage<double>(2); // two dimensional
    ar1->Allocate(numel,3);
    obj.GetEdips(di);
    for(int i1=0;i1<numel;i1++)
        {
            ar1->data2D[i1][0]=di->data1D[i1].x();
            ar1->data2D[i1][1]=di->data1D[i1].y();
            ar1->data2D[i1][2]=di->data1D[i1].z();
        }
    name = filenamed + ".edips";
    tio.WriteRectangular(ar1,name);
    ar1->Delete();
    delete ar1;
    di->Delete();
    delete di;


    // transition dipoles double
    storage<dvector3d>* dips2 = new storage<dvector3d>(2); // two dimensional
	dips2->Allocate(numel,numel2);
    storage<double>* ar = new storage<double>(2); // two dimensional
    ar->Allocate(numel*numel2,3);
    obj.GetDips2(dips2);
    for(int i1=0;i1<numel;i1++)
        for(int i2=0;i2<numel2;i2++)
        {
            ar->data2D[i1*numel2+i2][0]=dips2->data2D[i1][i2].x();
            ar->data2D[i1*numel2+i2][1]=dips2->data2D[i1][i2].y();
            ar->data2D[i1*numel2+i2][2]=dips2->data2D[i1][i2].z();
        }
    name = filenamed + ".dips2";
    tio.WriteRectangular(ar,name);
    
    obj.GetEdips2(dips2);
    for(int i1=0;i1<numel;i1++)
        for(int i2=0;i2<numel2;i2++)
        {
            ar->data2D[i1*numel2+i2][0]=dips2->data2D[i1][i2].x();
            ar->data2D[i1*numel2+i2][1]=dips2->data2D[i1][i2].y();
            ar->data2D[i1*numel2+i2][2]=dips2->data2D[i1][i2].z();
        }
    name = filenamed + ".edips2";
    tio.WriteRectangular(ar,name);
    ar->Delete();
    delete ar;
    dips2->Delete();
    delete dips2;
    // finally writing out pattern
    storage<double>* arp = new storage<double>(2); // two dimensional
    arp->Allocate(numel*numel2,2);
    for(int i1=0;i1<numel;i1++)
        for(int i2=0;i2<numel2;i2++)
        {
            arp->data2D[i1*numel2+i2][0]=i1;
            arp->data2D[i1*numel2+i2][1]=i2;
        }
    name = filenamed + ".dips2.pattern";
    tio.WriteRectangular(arp,name);
    arp->Delete();
    delete arp;
    
    
    
    

    // that is it.
}

