#include<iostream>

#include"dvectorNd.hpp"
#include<string>

using std::cout;

int main()
{
	dvectorNd vec1(10);
	dvectorNd vec2(10);

	vec2 = vec1;
	std::cout<<(vec1==vec2)<<"\n";

}
