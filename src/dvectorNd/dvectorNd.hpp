#pragma once

#include <sstream>
#include <string>


class dvectorNd
{
	public:

	//default constructor;
	dvectorNd();
	dvectorNd(int idim);
	dvectorNd(double* iv,int idim);
	//dvectorNd(double* iv,int idim,int);
	dvectorNd(const dvectorNd&);
	~dvectorNd();

	double Magnitude();
	double Amplitude();
	void Normalize();

	//operators
	dvectorNd& operator = (double rhs);// this is meaningful only when rhs = 0

	// assignment
	dvectorNd& operator = (const dvectorNd& rhs);
	dvectorNd& operator += (const dvectorNd& rhs);
	dvectorNd& operator -= (const dvectorNd& rhs);

	// assignment to value
	//dvectorNd& operator = (const double& rhs);
	//dvectorNd& operator = (const int& rhs);

	// multiplication by a scalar
	dvectorNd operator * (double rhs);
	//dvectorNd operator * (int rhs);

	// linear operations
	dvectorNd operator + (const dvectorNd& rhs);
	dvectorNd operator - (const dvectorNd& rhs);

	// scalar product
	double operator * (const dvectorNd& rhs);

	// comparison
	const bool operator == (const dvectorNd& rhs);
	const bool operator != (const dvectorNd& rhs);


	// communications
	//void Read(std::string);
	//std::string Write();

	double GetAt(int ind);

public:
	// do not use these in public unless necessary

	// the coordinates will be stored in
	//storage<double>* dcords;
	//void Allocate();
	//void Delete();

	int dim;

	// buffer is of dim size
	double* cooo;
	//bool bufferReady;

};

