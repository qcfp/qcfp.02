#pragma once


// include section
#include<iostream>
#include<cstring>
#include<typeinfo>



using std::cout;
using std::cout;

using std::cin;
using std::string;
using std::type_info;


//////////
// definition

// this is multidimensional array storage container
// it is designed so that the dimensionality should be defined and 
// it should not be changed. Otherwise too much confusion

// DO NOT DELETE ARRAYS OUTSIDE OF THE CLASS
// DO NOT CHANGE ELEMENTS OUTSIDE THE CLASS 
//     IMPLYING THAT THEY WILL BE SAVED

template <typename T>
class storage
{
public:
    // data storage
    // if you access them directly you have to know what you are doing
    // memory corruption can be easily done
    T* data1D;
    T** data2D;
    T*** data3D;
    T**** data4D;
    T***** data5D;
    T****** data6D;
    
private:
        // simple or complex type
    // this is needed when deciding how to do assignment
    int type_primitive;
    
    // dimension: 0 = unallocated
    int Dimension;
    
    // this can tell if the memory is alloated
    bool ready;
    
    // this can be used to check if elements can be changed
    bool readonly;
    
    // this can be used to check if elements can be accessed
    //bool locked;
    
    // number of elements in all dimensions
    // up to 6 dimensions
    int numEl[6];

    
    ////////////////
public:
    
    
	// creation with no dimension specification
	storage<T>()
    {
        data1D = 0;
        data2D = 0;
        data3D = 0;
        data4D = 0;
        data5D = 0;
        data6D = 0;
        
        Dimension = 0;
        
        numEl[0] = 0;
        numEl[1] = 0;
        numEl[2] = 0;
        numEl[3] = 0;
        numEl[4] = 0;
        numEl[5] = 0;
        
        ready = false;
        readonly = false;
        //locked = true;
        
        type_primitive = 0;
        if(typeid(T) == typeid(int))
            type_primitive = 1;
        if(typeid(T) == typeid(float))
            type_primitive = 1;
        if(typeid(T) == typeid(double))
            type_primitive = 1;
        
    }

	// creation with dimension but no allocation
	storage<T>(int dim)
    {
        data1D = 0;
        data2D = 0;
        data3D = 0;
        data4D = 0;
        data5D = 0;
        data6D = 0;
        
        Dimension = dim;
        
        numEl[0] = 0;
        numEl[1] = 0;
        numEl[2] = 0;
        numEl[3] = 0;
        numEl[4] = 0;
        numEl[5] = 0;
        
        ready = false;
        readonly = false;
        //locked = true;
        
        type_primitive = 0;
        if(typeid(T) == typeid(int))
            type_primitive = 1;
        if(typeid(T) == typeid(float))
            type_primitive = 1;
        if(typeid(T) == typeid(double))
            type_primitive = 1;
        
    }
	storage<T>(const storage<T>& rhs)
    {
        data1D = 0;
        data2D = 0;
        data3D = 0;
        data4D = 0;
        data5D = 0;
        data6D = 0;
        
        Dimension = 0;
        
        numEl[0] = 0;
        numEl[1] = 0;
        numEl[2] = 0;
        numEl[3] = 0;
        numEl[4] = 0;
        numEl[5] = 0;
        
        ready = false;
        readonly = false;
        //locked = true;
        
        type_primitive = 0;
        if(typeid(T) == typeid(int))
            type_primitive = 1;
        if(typeid(T) == typeid(float))
            type_primitive = 1;
        if(typeid(T) == typeid(double))
            type_primitive = 1;
        

        if(this == &rhs)
        {
            // skip
            return;
        }
        else if(&rhs == 0)
        {
            // skip
            return;
        }
        else if(!rhs.IsAlloc())
        {
            Dimension = rhs.Dimension;
        }
        else // normal copy
        {
            Allocate(rhs.Dimension,rhs.numEl);
            
            // assignment of data
            if(Dimension == 1)
                SetBar(rhs.data1D, rhs.numEl[0]);
            
            else if(Dimension==2)
            {
                for(int j1=0; j1<rhs.numEl[1]; j1++)
                    SetBar(j1, rhs.data2D[j1], rhs.numEl[0]);
            }
            
            else if(Dimension==3)
            {
                for(int j1=0; j1<rhs.numEl[2]; j1++)
                    for(int j2=0; j2<rhs.numEl[1]; j2++)
                        SetBar(j1, j2, rhs.data3D[j1][j2], rhs.numEl[0]);
            }
            else if(Dimension==4)
            {
                for(int j1=0; j1<rhs.numEl[3]; j1++)
                    for(int j2=0; j2<rhs.numEl[2]; j2++)
                        for(int j3=0; j3<rhs.numEl[1]; j3++)
                            SetBar(j1, j2, j3, rhs.data4D[j1][j2][j3], rhs.numEl[0]);
            }
            else if(Dimension==5)
            {
                for(int j1=0; j1<rhs.numEl[4]; j1++)
                    for(int j2=0; j2<rhs.numEl[3]; j2++)
                        for(int j3=0; j3<rhs.numEl[2]; j3++)
                            for(int j4=0; j4<rhs.numEl[1]; j4++)
                                SetBar(j1, j2, j3, j4, rhs.data5D[j1][j2][j3][j4], rhs.numEl[0]);
            }
            else if(Dimension==6)
            {
                for(int j1=0; j1<rhs.numEl[5]; j1++)
                    for(int j2=0; j2<rhs.numEl[4]; j2++)
                        for(int j3=0; j3<rhs.numEl[3]; j3++)
                            for(int j4=0; j4<rhs.numEl[2]; j4++)
                                for(int j5=0; j5<rhs.numEl[1]; j5++)
                                    SetBar(j1, j2, j3, j4, j5, rhs.data6D[j1][j2][j3][j4][j5], rhs.numEl[0]);
            }
            // assignment of flags
            type_primitive = rhs.type_primitive;
            
            ready = rhs.ready;
            
            readonly = rhs.readonly;
            
            //locked = rhs.locked;
        }
    }

	// destruction
	~storage<T>();



	// data creation routines
	void SetDimension( int dim)
	{ // only after default constructor it is posisble to call this function

		if(Dimension != 0)
			cout<<"Error: this function can be called only after default constructor!\n";
		
		if(dim == 0)
			cout<<"Error: construction with zero dimension is errorneous!\n";

		// otherwise it is OK
		Dimension = dim;
	}
	
	// 1D allocation
	int Allocate(int);
	// 2D allocation
	int Allocate(int,int);
	// 3D allocation
	int Allocate(int,int,int);
	// 4D allocation
	int Allocate(int,int,int,int);
	// 5D allocation
	int Allocate(int,int,int,int,int);
	// 6D allocation
	int Allocate(int,int,int,int,int,int);
	// deletion (un-allocation) of memory
	void Delete();

	// data access flag routined

	// this can tell if the memory is alloated
	bool IsAlloc() const {return ready;};
	// this can be used to check if elements can be changed
	bool IsRO() {return readonly;};
	// this can be used to check if elements can be accessed
	//bool IsLocked() {return locked;};

	// data saving and retrieving routines
	
	int CheckDimension(){return Dimension;};
	int  GetSize();
	void GetSize(int&);
	void GetSize(int&,int&);
	void GetSize(int&,int&,int&);
	void GetSize(int&,int&,int&,int&);
	void GetSize(int&,int&,int&,int&,int&);
	void GetSize(int&,int&,int&,int&,int&,int&);

	// saves onedimensional array of N elements
	void SetBar(T*,int N);
	void SetBar(int i2,T*,int N);
	void SetBar(int i3,int i2,T*,int N);
	void SetBar(int i4,int i3,int i2,T*,int N);
	void SetBar(int i5,int i4,int i3,int i2,T*,int N);
	void SetBar(int i6,int i5,int i4,int i3,int i2,T*,int N);

	// returns onedimensional array of int elements at specified place
	void GetBar(T*,int N);
	void GetBar(int i2,T*,int N);
	void GetBar(int i3,int i2,T*,int N);
	void GetBar(int i4,int i3,int i2,T*,int N);
	void GetBar(int i5,int i4,int i3,int i2,T*,int N);
	void GetBar(int i6,int i5,int i4,int i3,int i2,T*,int N);
	
	// for ultra fast access of the data
	// beware of possible memory corruption
	// know what you are doing: no checking is performed
	T* FlipBar(T*);
	T* FlipBar(int& i2,T*);
	T* FlipBar(int& i3,int& i2,T*);
	T* FlipBar(int& i4,int& i3,int& i2,T*);
	T* FlipBar(int& i5,int& i4,int& i3,int& i2,T*);
	T* FlipBar(int& i6,int& i5,int& i4,int& i3,int& i2,T*);
        
         // nD allocation
 	void Allocate(const int,const int*);
 	void GetSize(int*);
    
    const storage<T>& operator = (const storage<T>& rhs);
    const storage<T>& operator += (const storage<T>& rhs);

 
    void FillLinear(T valini,T valstep,int numpoints)
    {
        storage<T> tmps(1);
        tmps.Allocate(numpoints);
        
        for(int ind=0; ind<numpoints; ind++)
        {
            tmps.data1D[ind] = valini + valstep*ind;
        }
        *this = tmps;
    }
 	
 	// changes the size of the matrix
 	// exp is the expansion or substraction of the elements of the dimn dimension
 	// exp can be <0, dimn cannot
// 	void Resize(int dimn, int exp);

 
        T* Update1D(T* idat, int& num)
        {
                T* ret = data1D;
                data1D = idat;
	
                int iret = numEl[0];
                numEl[0] = num;
                num = iret;

                ready = true;
                readonly = false;
                //locked = false;
    
                return ret;
        }

private:
	T* Alloc(int);
	T** Alloc(int,int);
	T*** Alloc(int,int,int);
	T**** Alloc(int,int,int,int);
	T***** Alloc(int,int,int,int,int);
	T****** Alloc(int,int,int,int,int,int);
	void Del(T******,int,int,int,int,int,int);
	void Del(T*****,int,int,int,int,int);
	void Del(T****,int,int,int,int);
	void Del(T***,int,int,int);
	void Del(T**,int,int);
	void Del(T*,int);

	// the conditions to change data
	inline bool ModCChange()
	{
//		return (ready && !readonly && !locked );
        return (ready && !readonly );
	};
	// the conditions to read data
	inline bool ModCRead()
	{
//        return (ready && !locked );
		return (ready);
	};



};
////////
// implementation


////////////////////////////////////


template <typename T>
int storage<T>::Allocate(int n1)
{
    
    //cout<<Dimension<<" dimension\n";

    if(ready)
    {
        cout<<"Error: reallocation of allocated storage<T> is forbidden\n";
        return 1;
    }
	else if (Dimension == 1)
	{
            //cout<<n1<<" n1\n";
		data1D = Alloc(n1);
		numEl[0] = n1;
		ready = true;
		readonly = false;
		//locked = false;
        return 0;
	}
	else
		cout<<"Error: incompatible dimensions in allocation\n";

    return 2;
}
template <typename T>
int storage<T>::Allocate(int n2,int n1)
{
    if(ready)
    {
        cout<<"Error: reallocation of allocated storage<T> is forbidden\n";
        return 1;
    }
    else if (Dimension == 2)
	{
		data2D = Alloc(n2,n1);
		numEl[0] = n1;
		numEl[1] = n2;
		ready = true;
		readonly = false;
		//locked = false;
        return 0;
	}
	else
		cout<<"Error: incompatible dimensions in allocation\n";
    
    return 2;
}
template <typename T>
int storage<T>::Allocate(int n3,int n2,int n1)
{
    if(ready)
    {
        cout<<"Error: reallocation of allocated storage<T> is forbidden\n";
        return 1;
    }
    else if (Dimension == 3)
	{
		data3D = Alloc(n3,n2,n1);
		numEl[0] = n1;
		numEl[1] = n2;
		numEl[2] = n3;
		ready = true;
		readonly = false;
		//locked = false;
        return 0;
	}
	else
		cout<<"Error: incompatible dimensions in allocation\n";
    
    return 2;
}
template <typename T>
int storage<T>::Allocate(int n4,int n3,int n2,int n1)
{
    if(ready)
    {
        cout<<"Error: reallocation of allocated storage<T> is forbidden\n";
        return 1;
    }
    else if (Dimension == 4)
	{
		data4D = Alloc(n4,n3,n2,n1);
		numEl[0] = n1;
		numEl[1] = n2;
		numEl[2] = n3;
		numEl[3] = n4;
		ready = true;
		readonly = false;
        //locked = false;
        return 0;
	}
	else
		cout<<"Error: incompatible dimensions in allocation\n";
    
    return 2;
}
template <typename T>
int storage<T>::Allocate(int n5,int n4,int n3,int n2,int n1)
{
    if(ready)
    {
        cout<<"Error: reallocation of allocated storage<T> is forbidden\n";
        return 1;
    }
    else if (Dimension == 5)
	{
		data5D = Alloc(n5,n4,n3,n2,n1);
		numEl[0] = n1;
		numEl[1] = n2;
		numEl[2] = n3;
		numEl[3] = n4;
		numEl[4] = n5;
		ready = true;
		readonly = false;
		//locked = false;
        return 0;
	}
	else
		cout<<"Error: incompatible dimensions in allocation\n";
    
    return 2;
}
template <typename T>
int storage<T>::Allocate(int n6,int n5,int n4,int n3,int n2,int n1)
{
    if(ready)
    {
        cout<<"Error: reallocation of allocated storage<T> is forbidden\n";
        return 1;
    }
    else if (Dimension == 6)
	{
		data6D = Alloc(n6,n5,n4,n3,n2,n1);
		numEl[0] = n1;
		numEl[1] = n2;
		numEl[2] = n3;
		numEl[3] = n4;
		numEl[4] = n5;
		numEl[5] = n6;
		ready = true;
		readonly = false;
		//locked = false;
        return 0;
	}
	else
		cout<<"Error: incompatible dimensions in allocation\n";
    
    return 2;
}

template <typename T>
storage<T>::~storage()
{
	if(Dimension == 0)
		return;// do nothing
	
    if(ready)
        Delete();
		
    data1D = 0;
    data2D = 0;
    data3D = 0;
    data4D = 0;
    data5D = 0;
    data6D = 0;

    Dimension = 0;
    
    numEl[0] = 0;
    numEl[1] = 0;
    numEl[2] = 0;
    numEl[3] = 0;
    numEl[4] = 0;
    numEl[5] = 0;
        
    ready = false;
    readonly = false;
    //locked = true;
}


////////////////////////////////////
template <typename T>
void storage<T>::Delete()
{
    if(!ready)
        return;
    
    
	if (Dimension==0)
		return; // do nothing
	else if(Dimension == 1)
		Del(data1D,numEl[0]);
	else if(Dimension == 2)
		Del(data2D,numEl[1],numEl[0]);
	else if(Dimension == 3)
		Del(data3D,numEl[2],numEl[1],numEl[0]);
	else if(Dimension == 4)
		Del(data4D,numEl[3],numEl[2],numEl[1],numEl[0]);
	else if(Dimension == 5)
		Del(data5D,numEl[4],numEl[3],numEl[2],numEl[1],numEl[0]);
	else if(Dimension == 6)
		Del(data6D,numEl[5],numEl[4],numEl[3],numEl[2],numEl[1],numEl[0]);


    data1D = 0;
    data2D = 0;
    data3D = 0;
    data4D = 0;
    data5D = 0;
    data6D = 0;

    numEl[0] = 0;
    numEl[1] = 0;
    numEl[2] = 0;
    numEl[3] = 0;
    numEl[4] = 0;
    numEl[5] = 0;
		
    ready = false;
    readonly = false;
    //locked = true;

    // dimension remains the same
    // and the container can be reused
}


template <typename T>
void storage<T>::Del(T****** ar, int n6,int n5,int n4,int n3,int n2,int n1)
{
	for(int i=0;i<n6;i++)
		Del(ar[i],n5,n4,n3,n2,n1);
    if(n6 != 0)
	delete[] ar;
	ar = 0;
}
template <typename T>
void storage<T>::Del(T***** ar, int n5,int n4,int n3,int n2,int n1)
{
	for(int i=0;i<n5;i++)
		Del(ar[i],n4,n3,n2,n1);
    if(n5 != 0)
	delete[] ar;
	ar = 0;
}
template <typename T>
void storage<T>::Del(T**** ar, int n4,int n3,int n2,int n1)
{
	for(int i=0;i<n4;i++)
		Del(ar[i],n3,n2,n1);
    if(n4 != 0)
	delete[] ar;
	ar = 0;
}
template <typename T>
void storage<T>::Del(T*** ar, int n3,int n2,int n1)
{
	for(int i=0;i<n3;i++)
		Del(ar[i],n2,n1);
    if(n3 != 0)
	delete[] ar;
	ar = 0;
}
template <typename T>
void storage<T>::Del(T** ar, int n2,int n1)
{
	for(int i=0;i<n2;i++)
		Del(ar[i],n1);
    if(n2 != 0)
	delete[] ar;
	ar = 0;
}
template <typename T>
void storage<T>::Del(T* ar, int n1)
{
    if(n1 != 0)
        delete[] ar;
	ar = 0;
}



template <typename T>
T* storage<T>::Alloc(int n)
{
	T* a=new T[n];
	for(int ind=0;ind<n;ind++) 
		a[ind]=0;
	return a;
}

template <typename T>
T** storage<T>::Alloc(int n2,int n1)
{	
	T** a=new T*[n2];
	for(int ind=0;ind<n2;ind++)
 		a[ind]=Alloc(n1);
	return a;
}

template <typename T>
T*** storage<T>::Alloc(int n3,int n2,int n1)
{	
	T*** a=new T** [n3];
	for(int ind=0;ind<n3;ind++)
 		a[ind]=Alloc(n2,n1);
	return a;
}

template <typename T>
T**** storage<T>::Alloc(int n4,int n3,int n2,int n1)
{	
	T**** a=new T*** [n4];
	for(int ind=0;ind<n4;ind++)
 		a[ind]=Alloc(n3,n2,n1);
	return a;
}

template <typename T>
T***** storage<T>::Alloc(int n5,int n4,int n3,int n2,int n1)
{	
	T***** a=new T**** [n5];
	for(int ind=0;ind<n5;ind++)
 		a[ind]=Alloc(n4,n3,n2,n1);
	return a;
}

template <typename T>
T****** storage<T>::Alloc(int n6,int n5,int n4,int n3,int n2,int n1)
{	
	T****** a=new T***** [n6];
	for(int ind=0;ind<n6;ind++)
 		a[ind]=Alloc(n5,n4,n3,n2,n1);
	return a;
}

template <typename T>
int storage<T>::GetSize()
{
	if(Dimension==1)
		return numEl[0];
	else
		cout<<"Error: trying to access info at wrong dimension\n";
    return 0;
}
template <typename T>
void storage<T>::GetSize(int& num1)
{
	if(Dimension==1)
		num1 = numEl[0];
	else
		cout<<"Error: trying to access info at wrong dimension\n";
}
template <typename T>
void storage<T>::GetSize(int& num2,int& num1)
{
	if(Dimension==2){
		num1 = numEl[0];
		num2 = numEl[1];
	}else
		cout<<"Error: trying to access info at wrong dimension\n";
}
template <typename T>
void storage<T>::GetSize(int& num3,int& num2,int& num1)
{
	if(Dimension==3){
		num1 = numEl[0];
		num2 = numEl[1];
		num3 = numEl[2];
	}else
		cout<<"Error: trying to access info at wrong dimension\n";
}
template <typename T>
void storage<T>::GetSize(int& num4,int& num3,int& num2,int& num1)
{
	if(Dimension==4){
		num1 = numEl[0];
		num2 = numEl[1];
		num3 = numEl[2];
		num4 = numEl[3];
	}else
		cout<<"Error: trying to access info at wrong dimension\n";
}
template <typename T>
void storage<T>::GetSize(int& num5,int& num4,int& num3,int& num2,int& num1)
{
	if(Dimension==5){
		num1 = numEl[0];
		num2 = numEl[1];
		num3 = numEl[2];
		num4 = numEl[3];
		num5 = numEl[4];
	}else
		cout<<"Error: trying to access info at wrong dimension\n";
}
template <typename T>
void storage<T>::GetSize(int& num6,int& num5,int& num4,int& num3,int& num2,int& num1)
{
	if(Dimension==6){
		num1 = numEl[0];
		num2 = numEl[1];
		num3 = numEl[2];
		num4 = numEl[3];
		num5 = numEl[4];
		num6 = numEl[5];
	}else
		cout<<"Error: trying to access info at wrong dimension\n";
}

/////////////////////////

template <typename T>
void storage<T>::SetBar(T* dat,int N)
{
	if(Dimension!=1)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements:"<<N<<" "<<numEl[0]<<"\n";
	else if(ModCChange()){
            if(type_primitive)
                    memcpy(data1D,dat,N*sizeof(T));
            else
		for(int ind=0;ind<N;ind++)
			data1D[ind] = dat[ind];
	}
	else
		cout<<"Warning: cannot save data\n";
}
template <typename T>
void storage<T>::SetBar(int i2,T* dat,int N)
{
	if(Dimension!=2)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(ModCChange()){
            if(type_primitive)
                memcpy(data2D[i2],dat,N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			data2D[i2][i1] = dat[i1];
	}
	else
		cout<<"Warning: cannot save data\n";
}
template <typename T>
void storage<T>::SetBar(int i3,int i2,T* dat,int N)
{
	if(Dimension!=3)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(ModCChange()){
            if(type_primitive)
                memcpy(data3D[i3][i2],dat,N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			data3D[i3][i2][i1] = dat[i1];
	}
	else
		cout<<"Warning: cannot save data\n";
}
template <typename T>
void storage<T>::SetBar(int i4,int i3,int i2,T* dat,int N)
{
	if(Dimension!=4)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(i4<0 || i4>=numEl[3])
		cout<<"Error: wrong number of elements\n";
	else if(ModCChange()){
            if(type_primitive)
                memcpy(data4D[i4][i3][i2],dat,N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			data4D[i4][i3][i2][i1] = dat[i1];
	}
	else
		cout<<"Warning: cannot save data\n";
}
template <typename T>
void storage<T>::SetBar(int i5,int i4,int i3,int i2,T* dat,int N)
{
	if(Dimension!=5)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(i4<0 || i4>=numEl[3])
		cout<<"Error: wrong number of elements\n";
	else if(i5<0 || i5>=numEl[4])
		cout<<"Error: wrong number of elements\n";
	else if(ModCChange()){
            if(type_primitive)
                memcpy(data5D[i5][i4][i3][i2],dat,N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			data5D[i5][i4][i3][i2][i1] = dat[i1];
	}
	else
		cout<<"Warning: cannot save data\n";
}
template <typename T>
void storage<T>::SetBar(int i6,int i5,int i4,int i3,int i2,T* dat,int N)
{
	if(Dimension!=6)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(i4<0 || i4>=numEl[3])
		cout<<"Error: wrong number of elements\n";
	else if(i5<0 || i5>=numEl[4])
		cout<<"Error: wrong number of elements\n";
	else if(i6<0 || i6>=numEl[5])
		cout<<"Error: wrong number of elements\n";
	else if(ModCChange()){
            if(type_primitive)
                memcpy(data6D[i6][i5][i4][i3][i2],dat,N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			data6D[i6][i5][i4][i3][i2][i1] = dat[i1];
	}
	else
		cout<<"Warning: cannot save data\n";
}


/////////////////////
template <typename T>
void storage<T>::GetBar(T* dat,int N)
{
	if(Dimension!=1)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(ModCRead()){
            if(type_primitive)
                memcpy(dat,data1D,N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			dat[i1] = data1D[i1];
	}
	else
		cout<<"Warning: cannot read data\n";
}
template <typename T>
void storage<T>::GetBar(int i2,T* dat,int N)
{
	if(Dimension!=2)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(ModCRead()){
            if(type_primitive)
                memcpy(dat,data2D[i2],N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			dat[i1] = data2D[i2][i1];
	}
	else
		cout<<"Warning: cannot read data\n";
}
template <typename T>
void storage<T>::GetBar(int i3,int i2,T* dat,int N)
{
	if(Dimension!=3)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(ModCRead()){
            if(type_primitive)
                memcpy(dat,data3D[i3][i2],N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			dat[i1] = data3D[i3][i2][i1];
	}
	else
		cout<<"Warning: cannot read data\n";
}
template <typename T>
void storage<T>::GetBar(int i4,int i3,int i2,T* dat,int N)
{
	if(Dimension!=4)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(i4<0 || i4>=numEl[3])
		cout<<"Error: wrong number of elements\n";
	else if(ModCRead()){
            if(type_primitive)
                memcpy(dat,data4D[i4][i3][i2],N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			dat[i1] = data4D[i4][i3][i2][i1];
	}
	else
		cout<<"Warning: cannot read data\n";
}
template <typename T>
void storage<T>::GetBar(int i5,int i4,int i3,int i2,T* dat,int N)
{
	if(Dimension!=5)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(i4<0 || i4>=numEl[3])
		cout<<"Error: wrong number of elements\n";
	else if(i5<0 || i5>=numEl[4])
		cout<<"Error: wrong number of elements\n";
	else if(ModCRead()){
           if(type_primitive)
               memcpy(dat,data5D[i5][i4][i3][i2],N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			dat[i1] = data5D[i5][i4][i3][i2][i1];
	}
	else
		cout<<"Warning: cannot read data\n";
}
template <typename T>
void storage<T>::GetBar(int i6,int i5,int i4,int i3,int i2,T* dat,int N)
{
	if(Dimension!=6)
		cout<<"Error: wrong dimension\n";
	else if(N>numEl[0] || N<0)
		cout<<"Error: wrong number of elements\n";
	else if(i2<0 || i2>=numEl[1])
		cout<<"Error: wrong number of elements\n";
	else if(i3<0 || i3>=numEl[2])
		cout<<"Error: wrong number of elements\n";
	else if(i4<0 || i4>=numEl[3])
		cout<<"Error: wrong number of elements\n";
	else if(i5<0 || i5>=numEl[4])
		cout<<"Error: wrong number of elements\n";
	else if(i6<0 || i6>=numEl[5])
		cout<<"Error: wrong number of elements\n";
	else if(ModCRead()){
            if(type_primitive)
                memcpy(dat,data6D[i6][i5][i4][i3][i2],N*sizeof(T));
            else
		for(int i1=0;i1<N;i1++)
			dat[i1] = data6D[i6][i5][i4][i3][i2][i1];
	}
	else
		cout<<"Warning: cannot read data\n";
}


template <typename T>
T*  storage<T>::FlipBar(T* dat)
{
	T* dtemp = data1D;
	data1D = dat;
	return dtemp;
}
template <typename T>
T*  storage<T>::FlipBar(int& i2,T* dat)
{
	T* dtemp = data2D[i2];
	data2D[i2] = dat;
	return dtemp;
}
template <typename T>
T*  storage<T>::FlipBar(int& i3,int& i2,T* dat)
{
	T* dtemp = data3D[i3][i2];
	data3D[i3][i2] = dat;
	return dtemp;
}
template <typename T>
T*  storage<T>::FlipBar(int& i4,int& i3,int& i2,T* dat)
{
	T* dtemp = data4D[i4][i3][i2];
	data4D[i4][i3][i2] = dat;
	return dtemp;
}
template <typename T>
T*  storage<T>::FlipBar(int& i5,int& i4,int& i3,int& i2,T* dat)
{
	T* dtemp = data5D[i5][i4][i3][i2];
	data5D[i5][i4][i3][i2] = dat;
	return dtemp;
}
template <typename T>
T*  storage<T>::FlipBar(int& i6,int& i5,int& i4,int& i3,int& i2,T* dat)
{
	T* dtemp = data6D[i6][i5][i4][i3][i2];
	data6D[i6][i5][i4][i3][i2] =dat;
	return dtemp;
}


 template <typename T>
 void storage<T>::Allocate(const int n, const int* num)
 {
        Dimension = n;
 	if(n==1)		Allocate(num[0]);
 	else if(n==2)	Allocate(num[1], num[0]);
 	else if(n==3)	Allocate(num[2], num[1], num[0]);
 	else if(n==4)	Allocate(num[3], num[2], num[1], num[0]);
 	else if(n==5)	Allocate(num[4], num[3], num[2], num[1], num[0]);
 	else if(n==6)	Allocate(num[5], num[4], num[3], num[2], num[1], num[0]);
    else{
        cout<<"Error with storage<T> dimension\n";
        return;
    }
     
     ready = true;
     readonly = false;
     //locked = false;

        
 }

//////> 	double nu;????

 template <typename T>
 void storage<T>::GetSize(int* num)
 {
 	if(Dimension==1)		num[0]=GetSize();
 	else if(Dimension==2)	GetSize(num[1], num[0]);
 	else if(Dimension==3)	GetSize(num[2], num[1], num[0]);
 	else if(Dimension==4)	GetSize(num[3], num[2], num[1], num[0]);
 	else if(Dimension==5)	GetSize(num[4], num[3], num[2], num[1], num[0]);
 	else if(Dimension==6)	GetSize(num[5], num[4], num[3], num[2], num[1], num[0]);
 }
 
 template <typename T>
 const storage<T>& storage<T>:: operator = (const storage<T>& rhs)
 {
     if(&rhs == 0)
     {
         cout<<"Error: storage<T>:: operator = (const storage<T>& rhs): cannot assign zero-address object\n";
         *this =storage<T>();
     }
     if(this != &rhs)
     {
         // setting up all arrays with respect to dimensions
         //const int* tnum = rhs.numEl;
         
         if(Dimension == rhs.Dimension)
         {
             if(rhs.IsAlloc())
             {
                 if( memcmp ( numEl, rhs.numEl, Dimension*sizeof(int) ) == 0 )
                 {
                     // data storages are identical and data can be copied

                     // assignment of data
                     if(Dimension == 1)
                     {
                         SetBar(rhs.data1D, rhs.numEl[0]);
                     }
                     else if(Dimension==2)
                     {
                         for(int j1=0; j1<rhs.numEl[1]; j1++)
                             SetBar(j1, rhs.data2D[j1], rhs.numEl[0]);
                     }
                     
                     else if(Dimension==3)
                     {
                         for(int j1=0; j1<rhs.numEl[2]; j1++)
                             for(int j2=0; j2<rhs.numEl[1]; j2++)
                                 SetBar(j1, j2, rhs.data3D[j1][j2], rhs.numEl[0]);
                     }
                     else if(Dimension==4)
                     {
                         for(int j1=0; j1<rhs.numEl[3]; j1++)
                             for(int j2=0; j2<rhs.numEl[2]; j2++)
                                 for(int j3=0; j3<rhs.numEl[1]; j3++)
                                     SetBar(j1, j2, j3, rhs.data4D[j1][j2][j3], rhs.numEl[0]);
                     }
                     else if(Dimension==5)
                     {
                         for(int j1=0; j1<rhs.numEl[4]; j1++)
                             for(int j2=0; j2<rhs.numEl[3]; j2++)
                                 for(int j3=0; j3<rhs.numEl[2]; j3++)
                                     for(int j4=0; j4<rhs.numEl[1]; j4++)
                                         SetBar(j1, j2, j3, j4, rhs.data5D[j1][j2][j3][j4], rhs.numEl[0]);
                     }
                     else if(Dimension==6)
                     {
                         for(int j1=0; j1<rhs.numEl[5]; j1++)
                             for(int j2=0; j2<rhs.numEl[4]; j2++)
                                 for(int j3=0; j3<rhs.numEl[3]; j3++)
                                     for(int j4=0; j4<rhs.numEl[2]; j4++)
                                         for(int j5=0; j5<rhs.numEl[1]; j5++)
                                             SetBar(j1, j2, j3, j4, j5, rhs.data6D[j1][j2][j3][j4][j5], rhs.numEl[0]);
                     }
                     // assignment of flags
                     type_primitive = rhs.type_primitive;
                     ready = rhs.ready;
                     readonly = rhs.readonly;
                     //locked = rhs.locked;
                 }
                 else
                 {
                     // dimensions are identical but the number of elements is not identical
                     Delete();
                     Allocate(rhs.Dimension,rhs.numEl);
                     *this = rhs;
                 }
             }
             else
             {
                 // dimensions are the same,
                 // rhs is created but unallocated
                 Delete();
             }
         }
         else
         {
             Delete();

             // even dimensions are different
             if(rhs.IsAlloc())
             {
                 Allocate(rhs.Dimension,rhs.numEl);
                 *this = rhs;
             }
             else
             {
                 Dimension = rhs.Dimension;
             }
         }
     }
     
     return *this;
 }

 

template <typename T>
const storage<T>& storage<T>:: operator += (const storage<T>& rhs)
{
    //if(this != &rhs)
    //{
        const int* tnum = rhs.numEl;
        
        if(Dimension == rhs.Dimension && memcmp ( numEl, tnum, Dimension*sizeof(int) ) == 0 )
        {
            
            // assignment of data
            if(Dimension == 1)
                for(int ind = 0; ind<tnum[0];ind++)
                    this->data1D[ind]+=rhs.data1D[ind];
            
            else if(Dimension==2)
            {
                for(int j1=0; j1<tnum[1]; j1++)
                    for(int ind = 0; ind<tnum[0];ind++)
                        this->data2D[j1][ind]+=rhs.data2D[j1][ind];
            }
            
            else if(Dimension==3)
            {
                for(int j1=0; j1<tnum[2]; j1++)
                    for(int j2=0; j2<tnum[1]; j2++)
                        for(int ind = 0; ind<tnum[0];ind++)
                            this->data3D[j1][j2][ind]+=rhs.data3D[j1][j2][ind];
            }
            else if(Dimension==4)
            {
                for(int j1=0; j1<tnum[3]; j1++)
                    for(int j2=0; j2<tnum[2]; j2++)
                        for(int j3=0; j3<tnum[1]; j3++)
                            for(int ind = 0; ind<tnum[0];ind++)
                                this->data4D[j1][j2][j3][ind]+=rhs.data4D[j1][j2][j3][ind];
            }
            else if(Dimension==5)
            {
                for(int j1=0; j1<tnum[4]; j1++)
                    for(int j2=0; j2<tnum[3]; j2++)
                        for(int j3=0; j3<tnum[2]; j3++)
                            for(int j4=0; j4<tnum[1]; j4++)
                                for(int ind = 0; ind<tnum[0];ind++)
                                    this->data5D[j1][j2][j3][j4][ind]+=rhs.data5D[j1][j2][j3][j4][ind];
            }
            else if(Dimension==6)
            {
                for(int j1=0; j1<tnum[5]; j1++)
                    for(int j2=0; j2<tnum[4]; j2++)
                        for(int j3=0; j3<tnum[3]; j3++)
                            for(int j4=0; j4<tnum[2]; j4++)
                                for(int j5=0; j5<tnum[1]; j5++)
                                    for(int ind = 0; ind<tnum[0];ind++)
                                        this->data6D[j1][j2][j3][j4][j5][ind]+=rhs.data6D[j1][j2][j3][j4][j5][ind];
            }
            
            
        }
        else
        {
            cout<<"Error: increment in storage can be only of the same types and sizes\n";
        }
    //}
    return *this;
}


/*
 template <typename T>
 void storage<T>::Resize(int dimn, int exp)
 {
 	if(exp==0||Dimension==0||dimn==0)
 		return;
 	if(dimn>Dimension || dimn<0)
 		cout<<"Error: wrong dimension\n";
 	int *initial = new int [Dimension];
 	GetSize(initial);
 	storage<T> newT(Dimension);
 	newT.Allocate(Dimension, initial);
 	storage<T> temp(Dimension);
 	newT.Allocate(Dimension, initial);
 	if(Dimension==1)		{temp.data1D=newT.data1D; newT.data1D=data1D; data1D=temp.data1D;}
 	else if(Dimension==2)	{temp.data2D=newT.data2D; newT.data2D=data2D; data2D=temp.data2D;}
 	else if(Dimension==3)	{temp.data3D=newT.data3D; newT.data3D=data3D; data3D=temp.data3D;}
 	else if(Dimension==4)	{temp.data4D=newT.data4D; newT.data4D=data4D; data4D=temp.data4D;}
 	else if(Dimension==5)	{temp.data5D=newT.data5D; newT.data5D=data5D; data5D=temp.data5D;}
 	else if(Dimension==6)	{temp.data6D=newT.data6D; newT.data6D=data6D; data6D=temp.data6D;}
 	Delete();
 	//cout<<initial[dimn-1]<<"\t";
 	initial[dimn-1]+=exp;
 	//cout<<initial[dimn-1]<<"\n";
 	Allocate(newT.CheckDimension(), initial);
 	initial[dimn-1]-=exp;
 	if(Dimension==1)	data1D=newT.FlipBar(data1D);
 	if(Dimension==2)	for(int indb=0; ((dimn==2) ? indb<((initial[1]<initial[1]+exp) ? initial[1] : initial[1]+exp) : indb<initial[1]); indb++)
 						for(int inda=0; ((dimn==1) ? inda<((initial[0]<initial[0]+exp) ? initial[0] : initial[0]+exp) : inda<initial[0]); inda++)
 							data2D[indb][inda]=newT.data2D[indb][inda];
 	if(Dimension==3)	for(int indc=0; ((dimn==3) ? indc<((initial[2]<initial[2]+exp) ? initial[2] : initial[2]+exp) : indc<initial[2]); indc++)
 						for(int indb=0; ((dimn==2) ? indb<((initial[1]<initial[1]+exp) ? initial[1] : initial[1]+exp) : indb<initial[1]); indb++)
 						for(int inda=0; ((dimn==1) ? inda<((initial[0]<initial[0]+exp) ? initial[0] : initial[0]+exp) : inda<initial[0]); inda++)
 							data3D[indc][indb][inda]=newT.data3D[indc][indb][inda];
 	if(Dimension==4)	for(int indd=0; ((dimn==4) ? indd<((initial[3]<initial[3]+exp) ? initial[3] : initial[3]+exp) : indd<initial[3]); indd++)
 						for(int indc=0; ((dimn==3) ? indc<((initial[2]<initial[2]+exp) ? initial[2] : initial[2]+exp) : indc<initial[2]); indc++)
 						for(int indb=0; ((dimn==2) ? indb<((initial[1]<initial[1]+exp) ? initial[1] : initial[1]+exp) : indb<initial[1]); indb++)
 						for(int inda=0; ((dimn==1) ? inda<((initial[0]<initial[0]+exp) ? initial[0] : initial[0]+exp) : inda<initial[0]); inda++)
 							data4D[indd][indc][indb][inda]=newT.data4D[indd][indc][indb][inda];
 	if(Dimension==5)	for(int inde=0; ((dimn==5) ? inde<((initial[4]<initial[4]+exp) ? initial[4] : initial[4]+exp) : inde<initial[4]); inde++)
 						for(int indd=0; ((dimn==4) ? indd<((initial[3]<initial[3]+exp) ? initial[3] : initial[3]+exp) : indd<initial[3]); indd++)						for(int indc=0; ((dimn==3) ? indc<((initial[2]<initial[2]+exp) ? initial[2] : initial[2]+exp) : indc<initial[2]); indc++)
 						for(int indb=0; ((dimn==2) ? indb<((initial[1]<initial[1]+exp) ? initial[1] : initial[1]+exp) : indb<initial[1]); indb++)
 						for(int inda=0; ((dimn==1) ? inda<((initial[0]<initial[0]+exp) ? initial[0] : initial[0]+exp) : inda<initial[0]); inda++)
 							data5D[inde][indd][indc][indb][inda]=newT.data5D[inde][indd][indc][indb][inda];
 	if(Dimension==6)	for(int indf=0; ((dimn==6) ? indf<((initial[5]<initial[5]+exp) ? initial[5] : initial[5]+exp) : indf<initial[5]); indf++)
 						for(int inde=0; ((dimn==5) ? inde<((initial[4]<initial[4]+exp) ? initial[4] : initial[4]+exp) : inde<initial[4]); inde++)
 						for(int indd=0; ((dimn==4) ? indd<((initial[3]<initial[3]+exp) ? initial[3] : initial[3]+exp) : indd<initial[3]); indd++)
 						for(int indc=0; ((dimn==3) ? indc<((initial[2]<initial[2]+exp) ? initial[2] : initial[2]+exp) : indc<initial[2]); indc++)
 						for(int indb=0; ((dimn==2) ? indb<((initial[1]<initial[1]+exp) ? initial[1] : initial[1]+exp) : indb<initial[1]); indb++)
 						for(int inda=0; ((dimn==1) ? inda<((initial[0]<initial[0]+exp) ? initial[0] : initial[0]+exp) : inda<initial[0]); inda++)
 							data6D[indf][inde][indd][indc][indb][inda]=newT.data6D[indf][inde][indd][indc][indb][inda];
 }
*/
 
