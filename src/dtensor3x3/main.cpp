#include<iostream>

#include"dtensor3x3.hpp"

using std::cout;

int main()
{
	dtensor3x3 vec1;
	dtensor3x3 vec2;
	
	vec1 = dtensor3x3(1,1,1, 0,0,0,  2,2,2);
	vec2 = dtensor3x3(0,0,0, 1,1,1,  0,0,0);

	vec1 = vec1 + vec2;

	cout<< 	vec1 << "\n";
}
