#pragma once
#include"../dvectorNd/dvectorNd.hpp"
#include"../dvector3d/dvector3d.hpp"
#include<iostream>

class dtensor3x3 : public dvectorNd
{
	public:

	//default constructor;
	dtensor3x3();
	//~dvector3d();
	dtensor3x3(const dtensor3x3&);
	dtensor3x3(const dvector3d&, const dvector3d&);
	dtensor3x3(double ixx, double ixy, double ixz,
		double iyx, double iyy, double iyz,
		double izx, double izy, double izz);
    double xx();
    double xy();
    double xz();
    double yx();
    double yy();
    double yz();
    double zx();
    double zy();
    double zz();

	// assignment
	dtensor3x3& operator = (const double rhs);
	dtensor3x3& operator = (const dtensor3x3& rhs);
	dtensor3x3& operator += (const dtensor3x3& rhs);
	dtensor3x3& operator -= (const dtensor3x3& rhs);


	// multiplication by a scalar
	dtensor3x3 operator * ( double rhs);

	// linear operations
	dtensor3x3 operator + (const dtensor3x3& rhs);
	dtensor3x3 operator - (const dtensor3x3& rhs);


	friend std::ostream& operator<< (std::ostream& os,  dtensor3x3& inp);
};

std::ostream& operator<< (std::ostream& os,  dtensor3x3& inp);

