#pragma once
#include"../dvectorNd/dvectorNd.hpp"

class dvector3d: public dvectorNd
{
	public:

	//default constructor;
	dvector3d();
	//~dvector3d(){};
	dvector3d(const dvector3d&);
	dvector3d(double ix, double iy, double iz);
	dvector3d(int ix, int iy, int iz);

//	double Magnitude();
//	double Amplitude();
//	void Normalize();
    double x();
    double y();
    double z();
	void xyz(double& rx, double& ry, double& rz);

	// assignment
	dvector3d& operator = (const double rhs);
	dvector3d& operator = (const dvector3d& rhs);
	dvector3d& operator += (const dvector3d& rhs);
	dvector3d& operator -= (const dvector3d& rhs);


	// multiplication by a scalar
	dvector3d operator * ( double rhs);
	//dvector3d operator * ( int rhs);

	// linear operations
	dvector3d operator + (const dvector3d& rhs);
	dvector3d operator - (const dvector3d& rhs);

	// scalar product
	double operator * (const dvector3d& rhs);

	// vector product
	dvector3d operator ^ (const dvector3d& rhs);

	// comparison
	const bool operator == (const dvector3d& rhs);
	const bool operator != (const dvector3d& rhs);


};

