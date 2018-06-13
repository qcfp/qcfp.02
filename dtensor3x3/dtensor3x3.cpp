#include"dtensor3x3.hpp"
#include<cstring>

dtensor3x3::dtensor3x3()
:dvectorNd(9)
{
    //cooo[0] = 0;
    //cooo[1] = 0;
    //cooo[2] = 0;
}

dtensor3x3::dtensor3x3(double ixx, double ixy, double ixz,
double iyx, double iyy, double iyz,
double izx, double izy, double izz)
:dvectorNd(9)
{
    cooo[0] = ixx;
    cooo[1] = ixy;
    cooo[2] = ixz;

    cooo[3] = iyx;
    cooo[4] = iyy;
    cooo[5] = iyz;

    cooo[6] = izx;
    cooo[7] = izy;
    cooo[8] = izz;
}
dtensor3x3::dtensor3x3(const dtensor3x3& orig)
:dvectorNd(9)
{
	memcpy(cooo,orig.cooo,9*sizeof(double));
}

dtensor3x3::dtensor3x3(const dvector3d& o2, const dvector3d& o1)
:dvectorNd(9)
{
	for(int iu=0;iu<3;iu++)
	for(int iv=0;iv<3;iv++)
		cooo[iu*3+iv]=o2.cooo[iu]*o1.cooo[iv];
}


dtensor3x3& dtensor3x3::operator = (const double rhs)
{
	cooo[0]=rhs;
	cooo[1]=rhs;
	cooo[2]=rhs;
	cooo[3]=rhs;
	cooo[4]=rhs;
	cooo[5]=rhs;
	cooo[6]=rhs;
	cooo[7]=rhs;
	cooo[8]=rhs;
	return *this;    // Return ref for multiple assignment
}
dtensor3x3& dtensor3x3::operator = (const dtensor3x3& rhs)
{

	if (this != &rhs)
	{
        memcpy(cooo,rhs.cooo,9*sizeof(double));
	}
	return *this;    // Return ref for multiple assignment
}

dtensor3x3& dtensor3x3::operator += (const dtensor3x3& rhs)
{
	cooo[0]+=rhs.cooo[0];
	cooo[1]+=rhs.cooo[1];
	cooo[2]+=rhs.cooo[2];

	cooo[3]+=rhs.cooo[3];
	cooo[4]+=rhs.cooo[4];
	cooo[5]+=rhs.cooo[5];

	cooo[6]+=rhs.cooo[6];
	cooo[7]+=rhs.cooo[7];
	cooo[8]+=rhs.cooo[8];

	return *this;    // Return ref for multiple assignment
}

dtensor3x3& dtensor3x3::operator -= (const dtensor3x3& rhs)
{
	if (this != &rhs)
	{
	cooo[0]-=rhs.cooo[0];
	cooo[1]-=rhs.cooo[1];
	cooo[2]-=rhs.cooo[2];

	cooo[3]-=rhs.cooo[3];
	cooo[4]-=rhs.cooo[4];
	cooo[5]-=rhs.cooo[5];

	cooo[6]-=rhs.cooo[6];
	cooo[7]-=rhs.cooo[7];
	cooo[8]-=rhs.cooo[8];
        }
        else
        {
 		cooo[0]=cooo[1]=cooo[2]=cooo[3]=cooo[4]=cooo[5]=0;
 		cooo[6]=cooo[7]=cooo[8]=0;
	}
	return *this;    // Return ref for multiple assignment
}

dtensor3x3 dtensor3x3::operator+(const dtensor3x3 &rhs) 
{
	return dtensor3x3(*this) += rhs;
}
dtensor3x3 dtensor3x3::operator-(const dtensor3x3 &rhs)
{
	return dtensor3x3(*this) -= rhs;
}
dtensor3x3 dtensor3x3::operator*( double rhs)
{
	return dtensor3x3(cooo[0]*rhs, cooo[1]*rhs, cooo[2]*rhs,cooo[3]*rhs, cooo[4]*rhs, cooo[5]*rhs,cooo[6]*rhs, cooo[7]*rhs, cooo[8]*rhs);
}


double dtensor3x3::xx()
{
    return cooo[0];
}
double dtensor3x3::xy()
{
    return cooo[1];
}
double dtensor3x3::xz()
{
    return cooo[2];
}
double dtensor3x3::yx()
{
    return cooo[3];
}
double dtensor3x3::yy()
{
    return cooo[4];
}
double dtensor3x3::yz()
{
    return cooo[5];
}

double dtensor3x3::zx()
{
    return cooo[6];
}

double dtensor3x3::zy()
{
    return cooo[7];
}

double dtensor3x3::zz()
{
    return cooo[8];
}

std::ostream& operator<< (std::ostream& os,  dtensor3x3& inp)
{
	os<<"("<<inp.xx()<<" ";
	os<<inp.xy()<<" ";
	os<<inp.xz()<<" / ";
	os<<inp.yx()<<" ";
	os<<inp.yy()<<" ";
	os<<inp.yz()<<" / ";
	os<<inp.zx()<<" ";
	os<<inp.zy()<<" ";
	os<<inp.zz()<<")";

	return os;
}






