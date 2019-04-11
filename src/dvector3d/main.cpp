#include<iostream>

#include"dvector3d.hpp"

using std::cout;

int main()
{
	dvector3d vec1;
	dvector3d vec2;
	
	vec1 = dvector3d(1.0 ,1.0 ,1.0 );
	vec2 = dvector3d(2.0 ,2.0 ,2.0 );

	vec1 = vec1*3.0;

	dvector3d vec3 = vec1^vec2;

	double scpr = vec3*vec1;
	
	cout<< 	scpr <<" = " << 0 <<" ?\n";

}
