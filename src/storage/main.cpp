#include<iostream>

#include"storage.hpp"

int main()
{
	storage<int>* IntArray3D;
	
	IntArray3D = new storage<int>(3);
	
	IntArray3D->Allocate(3,3,3);

	// now ultrafast data manipulation

	// 1. I know that the lowest dimension has 3 elements
	int i3,i2,i1;
	int* bar = new int[3];
	bar[0]=1;
	bar[1]=2;
	bar[2]=3;

	// storing this bar in element [0][0][*]
	i3=0;i2=0;
	bar=IntArray3D->FlipBar(i3,i2,bar);

	// now bar is saved and replaced with previous content of storage

	// adding more information
	bar[0]=4;
	bar[1]=5;
	bar[2]=6;

	// storing this bar in element [2][0][*]
	i3=2;i2=0;
	bar=IntArray3D->FlipBar(i3,i2,bar);

	for(int ind3=0;ind3<3;ind3++)
	for(int ind2=0;ind2<3;ind2++)
	{
		// data copying
		IntArray3D->GetBar(ind3,ind2,bar,3);
		// printing out the content
		for(int ind1=0;ind1<3;ind1++)
			cout<<ind3<<" "<<ind2<<" "<<ind1<<" "<<bar[ind1]<<"\n";
	}


	//deleting bar
	delete[] bar;

	// deleting storage
	IntArray3D->Delete();

	// note that FlipBar never did data copying
	// only memory addresses have been replaced
	// also memory corruption did not occur
	// since bar had always different memory locations from storage elements
        
        // testing assignemnt
        storage<int> stint(1);
        storage<int> stins(1);
        stint.Allocate(3);
        stins.Allocate(3);
        stins.data1D[0]=1;
        stins.data1D[1]=2;
        stins.data1D[2]=3;
        
        stint = stins;
        
        // printing result
        cout<<stint.data1D[0]<<" "<<stint.data1D[1]<<" "<<stint.data1D[2]<<"\n";

        
 	
}
