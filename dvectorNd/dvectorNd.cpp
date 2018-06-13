#include"dvectorNd.hpp"
#include<cstring>
#include<cmath>
#include<iostream>

dvectorNd::dvectorNd()
{
	dim = 1;
	cooo = new double[1];
	cooo[0]=0.0;
}
dvectorNd::dvectorNd(int idim)
{
	if(idim>0 && idim < 1025)
	{
		dim = idim;
		cooo = new double[idim];
		memset(cooo,0,dim*sizeof(double));
	}
	else
	{
		dim = 0;
		cooo = 0;
		std::cout<<"Error: dvectorNd(int) size of dvector is specified incorrectly\n";
	}
}

dvectorNd::dvectorNd(double* iv,int idim)
{
	if(idim>0 && idim < 1025)
	{
		dim = idim;
		cooo = new double[idim];
		memcpy(cooo,iv,idim*sizeof(double));
	}
	else
	{
		dim = 0;
		cooo = 0;
		std::cout<<"Error: dvectorNd(int) size of dvector is specified incorrectly\n";
	}
}


dvectorNd::dvectorNd(const dvectorNd& idv)
{
	// create a copy
	dim = idv.dim;
	cooo = 0;
	if(idv.cooo !=0 )
	{
		cooo = new double[dim];
		memcpy(cooo,idv.cooo,dim*sizeof(double));
	}
}
dvectorNd::~dvectorNd()
{
	//Delete();
	if(cooo!=0)
		delete[] cooo;
	dim = 0;
	cooo = 0;
}

double dvectorNd::Magnitude()
{
	double ret = 0;
	for(int ind=0; ind<dim; ind++)
		ret = cooo[ind]*cooo[ind];
	return ret;
}
double dvectorNd::Amplitude()
{
	return sqrt(Magnitude());
}
void dvectorNd::Normalize()
{
	double length = Amplitude();
	if( length != 0.0 )
	{
		for(int ind=0; ind<dim; ind++)
			cooo[ind] /= length;
	}
}
dvectorNd& dvectorNd::operator = (const dvectorNd& rhs)
{
	if (this != &rhs)
	{
		if(dim != rhs.dim )
		{
			dim = rhs.dim;
			if(cooo != 0) delete[] cooo;
			if(rhs.cooo!=0) 
			{
				cooo = new double[dim];
			}
		}
		if(rhs.cooo!=0) 
		memcpy(cooo,rhs.cooo,dim*sizeof(double));
	}
	return *this;    // Return ref for multiple assignment
}

dvectorNd& dvectorNd::operator += (const dvectorNd& rhs)
{
	if (this != &rhs)
	{
		if(dim != rhs.dim )
		{
			std::cout<<"Error: vectors of different size cannot be combined\n";
			return *this;
		}
		// assignment		
		for(int ind=0; ind<dim; ind++)
			cooo[ind]+=rhs.cooo[ind];

	}
	return *this;    // Return ref for multiple assignment
}
dvectorNd& dvectorNd::operator -= (const dvectorNd& rhs)
{
	if (this != &rhs)
	{
		if(dim != rhs.dim )
		{
			std::cout<<"Error: vectors of different size cannot be combined\n";
			return *this;
		}
		
		// assignment		
		for(int ind=0; ind<dim; ind++)
			cooo[ind]-=rhs.cooo[ind];

	}
	return *this;    // Return ref for multiple assignment
}

dvectorNd dvectorNd::operator+(const dvectorNd &rhs) 
{
	return dvectorNd(*this) += rhs;
}
dvectorNd dvectorNd::operator-(const dvectorNd &rhs)
{
	return dvectorNd(*this) -= rhs;
}

// this is the dot product
double dvectorNd::operator*(const dvectorNd &rhs)
{

	if (this != &rhs)
	{
		if(dim != rhs.dim )
		{
			std::cout<<"Error: vectors of different size cannot be combined\n";
			return 0.0;
		}
		double result = 0;
		// assignment		
		for(int ind=0; ind<dim; ind++)
			result += cooo[ind]*rhs.cooo[ind];
		return result;
	}

	return Magnitude();
}


dvectorNd dvectorNd::operator*(double rhs)
{
	if(cooo!=0)
	{
		for(int ind=0; ind<dim; ind++)
			cooo[ind] *= rhs;
	}
	return *this;
}

const bool dvectorNd::operator==(const dvectorNd &rhs)
{
	if (this != &rhs)
	{
		if(dim != rhs.dim ) return false;
		if(rhs.cooo == 0 && cooo != 0) return false;
		if(cooo == 0 && rhs.cooo != 0) return false;
		if(rhs.cooo!=0)
		{
			// dimensions are equal and both nonzero
			for(int id=0;id<dim;id++)
				if(cooo[id]!=rhs.cooo[id]) return false;
		}
		// dimensions are equal and both zero
		// or both identical
	}

	return true;
}
const bool dvectorNd::operator!=(const dvectorNd &rhs)
{
	return !(*this==rhs);
}


double dvectorNd::GetAt(int ind)
{
		if(ind>=0 && ind<dim)
			return cooo[ind];
		
		std::cout<<"Error: request of an element out of bounds in dvectorNd\n";

		return 0.0;		
}


