#include"dvector3d.hpp"
#include"../storage/storage.hpp"
#include<cmath>

dvector3d::dvector3d()
:dvectorNd(3)
{
    //cooo[0] = 0;
    //cooo[1] = 0;
    //cooo[2] = 0;
}

dvector3d::dvector3d(double ix, double iy, double iz)
:dvectorNd(3)
{
//	Allocate();
//
//	double* dtemp = 0;
//	dtemp = cooo->FlipBar(dtemp);
//	dtemp[0]=ix;
//	dtemp[1]=iy;
//	dtemp[2]=iz;
//	dtemp = cooo->FlipBar(dtemp);
    cooo[0] = ix;
    cooo[1] = iy;
    cooo[2] = iz;

}
//dvector3d::dvector3d(int ix, int iy, int iz)
//{
//    cooo[0] = (double)ix;
//    cooo[1] = (double)iy;
//    cooo[2] = (double)iz;
//}
dvector3d::dvector3d(const dvector3d& orig)
:dvectorNd(3)
{
//	Allocate();
//
//    if(orig.cooo == 0)
//        return;
//    if(!orig.cooo->IsAlloc())
//        return;
//
//	*this = orig;
	memcpy(cooo,orig.cooo,3*sizeof(double));
//    cooo[0] = orig.cooo[0];
//    cooo[1] = orig.cooo[1];
//    cooo[2] = orig.cooo[2];
}

//void dvector3d::Allocate()
//{
//	cooo = new storage<double>(1);
//	cooo->Allocate(3);
//}
//void dvector3d::Delete()
//{
//	delete cooo;
//}



//double dvector3d::Magnitude()
//{
//	// getting coords
//	//double* dtemp = 0;
//	//dtemp = cooo->FlipBar(dtemp);
//	double ret = cooo[0]*cooo[0]+cooo[1]*cooo[1]+cooo[2]*cooo[2];
//	//dtemp = cooo->FlipBar(dtemp);
//	return ret;
//}
//double dvector3d::Amplitude()
//{
//	return sqrt(Magnitude());
//}
//void dvector3d::Normalize()
//{
//	// getting coords
//	//double* dtemp = 0;
//	//dtemp = cooo->FlipBar(dtemp);
//	double length = sqrt(cooo[0]*cooo[0]+cooo[1]*cooo[1]+cooo[2]*cooo[2]);
//	if( length != 0.0 )
//	{
//		cooo[0] /= length;
//		cooo[1] /= length;
//		cooo[2] /= length;
//	}
//	//dtemp = cooo->FlipBar(dtemp);
//}
dvector3d& dvector3d::operator = (const dvector3d& rhs)
{

	if (this != &rhs)
	{
        //if(rhs.cooo == 0)
        //    return *this;
        //if(!rhs.cooo->IsAlloc())
        //    return *this;
        memcpy(cooo,rhs.cooo,3*sizeof(double));
//        cooo[0]=rhs.cooo[0];
//        cooo[1]=rhs.cooo[1];
//		cooo[2]=rhs.cooo[2];
	}
	return *this;    // Return ref for multiple assignment
}
//dvector3d& dvector3d::operator = (int rhs)
//{ // this is meaningful only when rhs = 0
//	//if (this != &rhs)
//	{
//		// getting lhs data object
//		//double* ltemp = 0;
//		//ltemp = this->cooo->FlipBar(ltemp);
//        
//		cooo[0]=(double)rhs;
//		cooo[1]=(double)rhs;
//		cooo[2]=(double)rhs;
//        
//		// putting back
//		//ltemp = this->cooo->FlipBar(ltemp);
//        
//	}
//	return *this;    // Return ref for multiple assignment
//}

//dvector3d& dvector3d::operator = (double rhs)
//{// this is meaningful only when rhs = 0
//
//	//if (this != &rhs)
//	{
//		// getting lhs data object
//		//double* ltemp = 0;
//		//ltemp = this->cooo->FlipBar(ltemp);
//        
//		cooo[0]=rhs;
//		cooo[1]=rhs;
//		cooo[2]=rhs;
//        
//		// putting back
//		//ltemp = this->cooo->FlipBar(ltemp);
//        
//	}
//	return *this;    // Return ref for multiple assignment
//}


dvector3d& dvector3d::operator += (const dvector3d& rhs)
{
	if (this != &rhs)
	{
        //if(rhs.cooo == 0)
        //    return *this;
        //if(!rhs.cooo->IsAlloc())
        //    return *this;
        

        // getting lhs data object
		//double* ltemp = 0;
		//ltemp = this->cooo->FlipBar(ltemp);

		// getting rhs data object
		//double* rtemp = 0;
		//rtemp = rhs.cooo->FlipBar(rtemp);
		
		cooo[0]+=rhs.cooo[0];
		cooo[1]+=rhs.cooo[1];
		cooo[2]+=rhs.cooo[2];

		// putting back
		//ltemp = this->cooo->FlipBar(ltemp);
		//rtemp = rhs.cooo->FlipBar(rtemp);

	}
 	else
 	{
 		//double* temp = 0;
 		//temp = this->cooo->FlipBar(temp);
		
 		cooo[0]*=2;
 		cooo[1]*=2;
 		cooo[2]*=2;
 		
 		//temp = this->cooo->FlipBar(temp);
 	}
	return *this;    // Return ref for multiple assignment
}

dvector3d& dvector3d::operator -= (const dvector3d& rhs)
{
	if (this != &rhs)
	{

        //if(rhs.cooo == 0)
         //   return *this;
        //if(!rhs.cooo->IsAlloc())
        //    return *this;
        

        // getting lhs data object
	//	double* ltemp = 0;
	//	ltemp = this->cooo->FlipBar(ltemp);

		// getting rhs data object
	//	double* rtemp = 0;
	//	rtemp = rhs.cooo->FlipBar(rtemp);
		
		cooo[0]-=rhs.cooo[0];
		cooo[1]-=rhs.cooo[1];
		cooo[2]-=rhs.cooo[2];

		// putting back
		//ltemp = this->cooo->FlipBar(ltemp);
		//rtemp = rhs.cooo->FlipBar(rtemp);
        }
        else
        {
 		//double* temp = 0;
 		//temp = this->cooo->FlipBar(temp);
 		
 		cooo[0]=cooo[1]=cooo[2]=0;
 		
 		//temp = this->cooo->FlipBar(temp);
	}
	return *this;    // Return ref for multiple assignment
}

dvector3d dvector3d::operator+(const dvector3d &rhs) 
{
	return dvector3d(*this) += rhs;
}
dvector3d dvector3d::operator-(const dvector3d &rhs)
{
	return dvector3d(*this) -= rhs;
}
double dvector3d::operator*(const dvector3d &rhs)
{

	//if (this != &rhs)
	//{
        //if(rhs.cooo == 0)
        //    return 0;
        //if(!rhs.cooo->IsAlloc())
        //    return 0;
        

		// getting lhs data object
		//double* ltemp = 0;
		//ltemp = this->cooo->FlipBar(ltemp);

		// getting rhs data object
		//double* rtemp = 0;
		//rtemp = rhs.cooo->FlipBar(rtemp);

		double result = cooo[0]*rhs.cooo[0] + cooo[1]*rhs.cooo[1] + cooo[2]*rhs.cooo[2];

		// putting back
		//ltemp = this->cooo->FlipBar(ltemp);
		//rtemp = rhs.cooo->FlipBar(rtemp);

		return result;
	//}

	// getting lhs data object
	//double* ltemp = 0;
	//ltemp = this->cooo->FlipBar(ltemp);

	//double result = ltemp[0]*ltemp[0] + ltemp[1]*ltemp[1] + ltemp[2]*ltemp[2];

	// putting back
	//ltemp = this->cooo->FlipBar(ltemp);

	//return result;
}

dvector3d dvector3d::operator^(const dvector3d &rhs) 
{

	//if (this != &rhs)
	//{

        //if(rhs.cooo == 0)
        //    return *this;
        //if(!rhs.cooo->IsAlloc())
        //    return *this;
        

        // getting lhs data object
		//double* ltemp = cooo;
		//ltemp = this->cooo->FlipBar(ltemp);

		// getting rhs data object
		//double* rtemp = rhs.cooo;
		//rtemp = rhs.cooo->FlipBar(rtemp);

		dvector3d temp( cooo[1]*rhs.cooo[2]-rhs.cooo[1]*cooo[2], cooo[2]*rhs.cooo[0]-rhs.cooo[2]*cooo[0], cooo[0]*rhs.cooo[1]-rhs.cooo[0]*cooo[1] );

		// putting back
		//ltemp = this->cooo->FlipBar(ltemp);
		//rtemp = rhs.cooo->FlipBar(rtemp);

		return temp;
	//}

	//return dvector3d();
}
dvector3d dvector3d::operator*( double rhs)
{
	return dvector3d(cooo[0]*rhs, cooo[1]*rhs, cooo[2]*rhs);
}
//
//dvector3d dvector3d::operator*( int rhs) 
//{
//	// getting lhs data object
//	//double* ltemp = 0;
//	//ltemp = this->cooo->FlipBar(ltemp);
//
//	dvector3d temp(cooo[0]*rhs, cooo[1]*rhs, cooo[2]*rhs);
//
//	// putting back
//	//ltemp = this->cooo->FlipBar(ltemp);
//
//	return temp;
//}

dvector3d& dvector3d::operator = (const double rhs)
{
	cooo[0]=rhs;
	cooo[1]=rhs;
	cooo[2]=rhs;
	return *this;
}

const bool dvector3d::operator==(const dvector3d &rhs)
{
	if (this != &rhs)
	{

		bool temp = ((cooo[0]==rhs.cooo[0]) && (cooo[1]==rhs.cooo[1]) && (cooo[2]==rhs.cooo[2]));

		return temp;
	}

	return true;
}
const bool dvector3d::operator!=(const dvector3d &rhs)
{
	return !(*this==rhs);
}


double dvector3d::x()
{
    return cooo[0];
}

double dvector3d::y()
{
    return cooo[1];
}

double dvector3d::z()
{
    return cooo[2];
}

void dvector3d::xyz(double& rx, double& ry, double& rz)
{
    rx = cooo[0];
    ry = cooo[1];
    rz = cooo[2];
}


