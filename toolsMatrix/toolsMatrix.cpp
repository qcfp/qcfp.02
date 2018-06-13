
#include"../storage/storage.hpp"
#include"../toolsMatrix/toolsMatrix.hpp"
#include<complex>
#include<cmath>
#include <stdlib.h>     /* malloc, free, rand */

#include"include-LAPACK/zgetrf.h"
#include"include-LAPACK/zgetri.h"


toolsMatrix::toolsMatrix()
{
    trd = 0;
    trdallocated = 0;

}
toolsMatrix::toolsMatrix(toolsRandom* itrd)
{
    trd = itrd;
    trdallocated = 2;
    
}
toolsMatrix::toolsMatrix(int seed)
{
    trd = new toolsRandom(seed);
    trdallocated = 1;
}
toolsMatrix::~toolsMatrix()
{
    if(trdallocated == 1) delete trd;
}




void toolsMatrix::Set0(storage<double>& matr)
{
  int il,ir;
  matr.GetSize(il,ir);
  for(int iil = 0; iil< il; iil++)
  for(int iir = 0; iir< ir; iir++)
      matr.data2D[iil][iir] = 0.0;
    
}


void toolsMatrix::Set1(storage<double>& matr)
{
  int il,ir;
  matr.GetSize(il,ir);
  for(int iil = 0; iil< il; iil++)
  for(int iir = 0; iir< ir; iir++)
  {
      if (iil == iir)
      matr.data2D[iil][iir] = 1.0;
      else
      matr.data2D[iil][iir] = 0.0;
  }
}





void toolsMatrix::VecRandLD(storage<double>& vec, double k)
//sets elements linear random from 0 to <k
{
    if(trdallocated)
    {
        int il;
        vec.GetSize(il);
        for(int iil = 0; iil< il; iil++)
        {
            vec.data1D[iil] = trd->RandLD(k);
        }
    }
    else
    {
        std::cout<<"Error: insufficient allocation of toolsMatrix for VecRandLD()\n";
    }
    
}
void toolsMatrix::VecRandED(storage<double>& vec, double k)
//sets elements exponential random
{
    if(trdallocated)
    {
        int il;
        vec.GetSize(il);
        for(int iil = 0; iil< il; iil++)
        {
            vec.data1D[iil] = trd->RandED(k);
        }
    }
    else
    {
        std::cout<<"Error: insufficient allocation of toolsMatrix for VecRandLD()\n";
    }
    
}
void toolsMatrix::VecRandGD(storage<double>& vec, double v, double s)
//sets elements gaussian random
{
    if(trdallocated)
    {
        int il;
        vec.GetSize(il);
        for(int iil = 0; iil< il; iil++)
        {
            vec.data1D[iil] = trd->RandGD(v,s);
        }
        
    }
    else
    {
        std::cout<<"Error: insufficient allocation of toolsMatrix for VecRandLD()\n";
    }
}









void toolsMatrix::AddVecToDiag(storage<double>& matr,storage<double>& vec)
{
    // adding a vector to diagonal
    // the length of a vector must be equal to the side-length of the matrix
    int il,ir;
    matr.GetSize(il,ir);
    
    for(int iil = 0; iil< ir; iil++)
    {
          matr.data2D[iil][iil] += vec.data1D[iil];
    }
    
}

void toolsMatrix::Prod (storage<double>& ret,storage<double>& mal,storage<double>& mar)
{
// this is the matrix multiplication
    // user must maintain proper matrix sizing
      int il,ir;
      mar.GetSize(il,ir);

  for(int iil = 0; iil< ir; iil++)
  for(int iim = 0; iim< il; iim++)
  for(int iir = 0; iir< ir; iir++)
  {
      ret.data2D[iil][iir] = mal.data2D[iil][iim]*mar.data2D[iim][iir];
  }
      
}

void toolsMatrix::Add (storage<double>& ret, storage<double>& mal, storage<double>& mar)
{
// this is the matrix addition
    // user must maintain proper matrix sizing
      int il,ir;
      mar.GetSize(il,ir);

  for(int iil = 0; iil< il; iil++)
  for(int iir = 0; iir< ir; iir++)
  {
      ret.data2D[iil][iir] = mal.data2D[iil][iir]+mar.data2D[iil][iir];
  }
      
}
void toolsMatrix::InverseC (storage<complexv>& ret, storage<complexv>& mat)
{
    // this is the matrix inversion
    int NN;
    mat.GetSize(NN);
    
    std::complex<double>* matrixIO = new std::complex<double>[NN*NN];
    for(int il = 0; il<NN; il++)
    for(int ir = 0; ir<NN; ir++)
        matrixIO[il*NN+ir] = mat.data2D[il][ir];
    
		
	//inverting matrix using LAPACK:
	int lapack_info = 0;
	int* lapack_iPivot  = (int*)malloc(NN*sizeof(int));
	std::complex<double>* lapack_cwork  = (std::complex<double>*)malloc(64*2*NN*sizeof(std::complex<double>));

//#ifdef __GNUC__
	zgetrf_(NN, NN, matrixIO, NN, lapack_iPivot, lapack_info);
	zgetri_(NN, matrixIO, NN, lapack_iPivot, lapack_cwork, 2*64*NN, lapack_info);

	free(lapack_iPivot);//  = (int*)malloc(NN*sizeof(int));
	free(lapack_cwork);//  = (complex<double>*)malloc(64*2*NN*sizeof(complex<double>));

    for(int il = 0; il<NN; il++)
    for(int ir = 0; ir<NN; ir++)
         ret.data2D[il][ir] = matrixIO[il*NN+ir];
    
    delete[] matrixIO;
}


void toolsMatrix::InverseC (complexv* ret, complexv* mat, int n)
{
    int NN=n;
    
    std::complex<double>* matrixIO = new std::complex<double>[NN*NN];
    for(int il = 0; il<NN*NN; il++)
        matrixIO[il] = mat[il];
		
	//inverting matrix using LAPACK:
	int lapack_info = 0;
	int* lapack_iPivot  = (int*)malloc(NN*sizeof(int));
	std::complex<double>* lapack_cwork  = (std::complex<double>*)malloc(64*2*NN*sizeof(std::complex<double>));

//#ifdef __GNUC__
	zgetrf_(NN, NN, matrixIO, NN, lapack_iPivot, lapack_info);
	zgetri_(NN, matrixIO, NN, lapack_iPivot, lapack_cwork, 2*64*NN, lapack_info);

	free(lapack_iPivot);//  = (int*)malloc(NN*sizeof(int));
	free(lapack_cwork);//  = (complex<double>*)malloc(64*2*NN*sizeof(complex<double>));

        for(int ii = 0; ii<NN*NN; ii++)
                mat[ii] = matrixIO[ii];
        if(ret != 0)
                for(int ii = 0; ii<NN*NN; ii++)
                        ret[ii] = matrixIO[ii];
    
        delete[] matrixIO;

}


