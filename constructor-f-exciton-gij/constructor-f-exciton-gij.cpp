// everything is performed
#include"constructor-f-exciton-gij.hpp"


#include"../eigen/eigen.hpp"
#include"../constants/constants.hpp"


constructor_f_exciton_gij::constructor_f_exciton_gij(storage<double>* ievecs)
{
    evecs = new storage<double>(2);
    
    *evecs = *ievecs;
    
//    int num;
//    evecs->GetSize(num,num);
//    for(int i2=0; i2<num; i2++)
//    for(int i1=0; i1<num; i1++)
//    {
//        cout<<i2<<" "<<i1<<" "<<evecs->data2D[i2][i1]<<" "<<ievecs->data2D[i2][i1]<<"\n";
//    }
    
    evec2 = 0;

}

constructor_f_exciton_gij::constructor_f_exciton_gij(storage<double>& ievecs)
{
    evecs = new storage<double>(2);
    *evecs = ievecs;
    evec2 = 0;
}

constructor_f_exciton_gij::constructor_f_exciton_gij(storage<double>* ievecs,storage<double>* ievec2)
{
    evecs = new storage<double>(2);
    *evecs = *ievecs;
    evec2 = 0;
   
    if(ievec2 !=0)
    {
	    evec2 = new storage<double>(3);
    	*evec2 = *ievec2;
        // modifying to make easier programming
        int numf,nume;
        evec2->GetSize(numf,nume,nume);
        for(int inf=0;inf<numf;inf++)
        for(int ind=0;ind<nume;ind++)
        {
            evec2->data3D[inf][ind][ind] *= constants::root2;
        }

    }
}

constructor_f_exciton_gij::constructor_f_exciton_gij(storage<double>& ievecs,storage<double>& ievec2)
{
    evecs = new storage<double>(2);
    *evecs = ievecs;

    evec2 = new storage<double>(3);
    *evec2 = ievec2;
    
    if(ievec2.IsAlloc())
    {
      // modifying to make easier programming
      int numf,nume;
      evec2->GetSize(numf,nume,nume);
      for(int inf=0;inf<numf;inf++)
      for(int ind=0;ind<nume;ind++)
            evec2->data3D[inf][ind][ind] *= constants::root2;
    }

}

constructor_f_exciton_gij::~constructor_f_exciton_gij()
{
    // no deallocation is done
    delete evecs;
    if(evec2 != 0)
	    delete evec2;
}

// cstrength is a number
int constructor_f_exciton_gij::GetGijAmplitudes11(storage<double>* ret,double& cstrength)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ia][sn]
                            *lev[ib][sn]
                            *lev[ib][sn];
            }
            amp[ia][ib] = cstrength*dtemp;
        }
    return 0;

}

// cstrength is a 1D vector
int constructor_f_exciton_gij::GetGijAmplitudes11(storage<double>* ret,double* cstrength)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ia][sn]
                            *lev[ib][sn]
                            *lev[ib][sn]
                        *cstrength[sn];
            }
            amp[ia][ib] = dtemp;
        }
    return 0;
}
    

//------------------------

// cstrength is a number
int constructor_f_exciton_gij::GetRedAmplitudes11(storage<double>* ret,double& cstrength)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ib][sn]
                            *lev[ib][sn]
                            *lev[ia][sn];
            }
            amp[ia][ib] = cstrength*dtemp;
        }

    return 0;
}

// cstrength is a 1D vector
int constructor_f_exciton_gij::GetRedAmplitudes11(storage<double>* ret,double* cstrength)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ib][sn]
                            *lev[ib][sn]
                            *lev[ia][sn]
                        *cstrength[sn];
            }
            amp[ia][ib] = dtemp;
        }
    return 0;
}
    
// cstrength is a 1D vector
int constructor_f_exciton_gij::GetGijAmplitudes21(storage<double>* ret,double* cstrength)
{
    double** lev = evecs->data2D;
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ie=0; ie<num; ie++ )
    {
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
            dtemp +=      lev[ie][ia]
                            *lev[ie][ia]
                            *le2[is][ia][ib]
                            *le2[is][ia][ib]
                        *cstrength[ia];
        }
        amp[is][ie] = dtemp;
    }
    return 0;

}
    
int constructor_f_exciton_gij::GetGijAmplitudes21(storage<double>* ret,double& cstrength)
{
    double** lev = evecs->data2D;
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ie=0; ie<num; ie++ )
    {
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
//            if(ia == ib )
//               dtemp +=    2*lev[ie][ia]
//                            *lev[ie][ia]
//                            *le2[is][ia][ib]
//                            *le2[is][ia][ib];
//            else
               dtemp +=      lev[ie][ia]
                            *lev[ie][ia]
                            *le2[is][ia][ib]
                            *le2[is][ia][ib];
        }
        amp[is][ie] = dtemp*cstrength;
    }
    return 0;

}
    
int constructor_f_exciton_gij::GetGijAmplitudes22(storage<double>* ret,double* cstrength)
{
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    { 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        {
//            if(ia != ib && ia != ic)
               dtemp +=      le2[is][ia][ib]
                            *le2[is][ia][ib]
                            *le2[ih][ia][ic]
                            *le2[ih][ia][ic]
                        *cstrength[ia];
//            else if(ia != ib && ia == ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[is][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[ih][ia][ic]
//                        *cstrength[ia];
//            else if(ia == ib && ia != ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[is][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[ih][ia][ic]
//                        *cstrength[ia];
//            else
 //              dtemp +=    4*le2[is][ia][ib]
 //                           *le2[is][ia][ib]
 //                           *le2[ih][ia][ic]
 //                           *le2[ih][ia][ic]
 //                       *cstrength[ia];
        }
        amp[is][ih] = dtemp;
    }

    return 0;
}
    
int constructor_f_exciton_gij::GetGijAmplitudes22(storage<double>* ret,double& cstrength)
{
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    { 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        {
//            if(ia != ib && ia != ic)
               dtemp +=      le2[is][ia][ib]
                            *le2[is][ia][ib]
                            *le2[ih][ia][ic]
                            *le2[ih][ia][ic];
//            else if(ia != ib && ia == ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[is][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[ih][ia][ic];
//            else if(ia == ib && ia != ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[is][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[ih][ia][ic];
//            else
//               dtemp +=    4*le2[is][ia][ib]
//                            *le2[is][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[ih][ia][ic];
        }
        amp[is][ih] = dtemp*cstrength;
    }
    return 0;

}
    
int constructor_f_exciton_gij::GetRedAmplitudes22(storage<double>* ret,double& cstrength)
{
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    { 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        {
//            if(ia != ib && ia != ic)
               dtemp +=      le2[is][ia][ib]
                            *le2[ih][ia][ib]
                            *le2[ih][ia][ic]
                            *le2[is][ia][ic];
//            else if(ia != ib && ia == ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[is][ia][ic];
//            else if(ia == ib && ia != ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[is][ia][ic];
//            else
//               dtemp +=    4*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[is][ia][ic];
        }
        amp[is][ih] = dtemp*cstrength;
    }
    return 0;

}
    
int constructor_f_exciton_gij::GetRedAmplitudes22(storage<double>* ret,double* cstrength)
{
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    { 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        {
//            if(ia != ib && ia != ic)
               dtemp +=      le2[is][ia][ib]
                            *le2[ih][ia][ib]
                            *le2[ih][ia][ic]
                            *le2[is][ia][ic]
                        *cstrength[ia];
//            else if(ia != ib && ia == ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[is][ia][ic]
//                        *cstrength[ia];
//            else if(ia == ib && ia != ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[is][ia][ic]
//                        *cstrength[ia];
//            else
//               dtemp +=    4*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[ih][ia][ic]
//                            *le2[is][ia][ic]
//                        *cstrength[ia];
        }
        amp[is][ih] = dtemp;
    }
    return 0;

}

int constructor_f_exciton_gij::GetFullRedAmplitudes11(storage<double>& R, storage<double>& cst)
{
    if(cst.CheckDimension()==1)
        return GetFullRedAmplitudes11(R.data4D, cst.data1D);
    else
        return GetFullRedAmplitudes11(R.data4D, cst.data2D);
}
int constructor_f_exciton_gij::GetFullRedAmplitudes11(double**** amp, double* cstrength)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        for( int id=0; id<num; id++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ib][sn]
                            *lev[ic][sn]
                            *lev[id][sn]
                        *cstrength[sn];
            }
            amp[ia][ib][ic][id] = dtemp;
        }
    return 0;
}
int constructor_f_exciton_gij::GetFullRedAmplitudes22(storage<double>& R, storage<double>& cst)
{
    if(cst.CheckDimension()==1)
        return GetFullRedAmplitudes22(R.data4D, cst.data1D);
    else
        return GetFullRedAmplitudes22(R.data4D, cst.data2D);
} 
int constructor_f_exciton_gij::GetFullRedAmplitudes22(double**** amp, double* cstrength)
{
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    for( int id=0; id<nu2; id++ )
    for( int ij=0; ij<nu2; ij++ )
    { 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        {
//            if(ia != ib && ia != ic)
               dtemp +=      le2[is][ia][ib]
                            *le2[ih][ia][ib]
                            *le2[id][ia][ic]
                            *le2[ij][ia][ic]
                        *cstrength[ia];
//            else if(ia != ib && ia == ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[id][ia][ic]
//                            *le2[ij][ia][ic]
//                        *cstrength[ia];
//            else if(ia == ib && ia != ic)
//               dtemp +=    2*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[id][ia][ic]
//                            *le2[ij][ia][ic]
//                        *cstrength[ia];
//            else
//               dtemp +=    4*le2[is][ia][ib]
//                            *le2[ih][ia][ib]
//                            *le2[id][ia][ic]
//                            *le2[ij][ia][ic]
//                        *cstrength[ia];
        }
        amp[is][ih][id][ij] = dtemp;
    }
    return 0;
}

int constructor_f_exciton_gij::GetFullRedAmplitudes21(storage<double>& R, storage<double>& cst)
{
    if(cst.CheckDimension()==1)
        return GetFullRedAmplitudes21(R.data4D, cst.data1D);
    else
        return GetFullRedAmplitudes21(R.data4D, cst.data2D);
}
int constructor_f_exciton_gij::GetFullRedAmplitudes21(double**** amp, double* cstrength)
{
    double** lev = evecs->data2D;
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int iz=0; iz<nu2; iz++ )
    for( int ie=0; ie<num; ie++ )
    for( int ix=0; ix<num; ix++ )
    {
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
               dtemp +=      lev[ie][ia]
                            *lev[ix][ia]
                            *le2[is][ia][ib]
                            *le2[iz][ia][ib]
                            *cstrength[ia];
        }
        amp[is][iz][ie][ix] = dtemp;
    }
    return 0;
}



//////////////////////////
int constructor_f_exciton_gij::GetGijAmplitudesC11(storage<double>& ret, double** cst)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;

    double** amp = ret.data2D;
    double** ccc = cst;
    
    // here we assume that all sites fluctuate with the same type of fluctuation with correlations.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            for( int sm=0; sm<num; sm++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ia][sn]
                            *lev[ib][sm]
                            *lev[ib][sm]
                        *ccc[sn][sm];
            }
            amp[ia][ib] = dtemp;
        }
    return 0;
}

int constructor_f_exciton_gij::GetRedAmplitudesC11(storage<double>& ret, double** cst)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;

    double** amp = ret.data2D;
    double** ccc = cst;
    
    // here we assume that all sites fluctuate with the same type of fluctuation with correlations.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            for( int sm=0; sm<num; sm++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ib][sn]
                            *lev[ib][sm]
                            *lev[ia][sm]
                        *ccc[sn][sm];
            }
            amp[ia][ib] = dtemp;
        }
    return 0;
}
    
int constructor_f_exciton_gij::GetGijAmplitudesC21(storage<double>& ret, double** cst)
{
    double** lev = evecs->data2D;
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret.data2D;
    double** ccc = cst;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ie=0; ie<num; ie++ )
    {
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        {
//            if(ia == ib )
//               dtemp +=    2*lev[ie][ia]
//                            *lev[ie][ia]
//                            *le2[is][ia][ib]
//                            *le2[is][ia][ib];
//            else
               dtemp +=      lev[ie][ia]
                            *lev[ie][ia]
                            *le2[is][ic][ib]
                            *le2[is][ic][ib]
                        *ccc[ia][ib];
        }
        amp[is][ie] = dtemp;
    }
    return 0;

}

int constructor_f_exciton_gij::GetGijAmplitudesC22(storage<double>& ret, double** cst)
{
    double** lev = evecs->data2D;
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret.data2D;
    double** ccc = cst;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
   for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    { 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        for( int id=0; id<num; id++ )
        {
               dtemp +=      le2[is][ia][ib]
                            *le2[is][ia][ib]
                            *le2[ih][id][ic]
                            *le2[ih][id][ic]
                        *ccc[ia][ic];
        }
        amp[is][ih] = dtemp;
    }
    return 0;

}

int constructor_f_exciton_gij::GetRedAmplitudesC22(storage<double>& ret, double** cst)
{
    double** lev = evecs->data2D;
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    double** amp = ret.data2D;
    double** ccc = cst;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
   for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    { 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        for( int id=0; id<num; id++ )
        {
               dtemp +=      le2[is][ia][ib]
                            *le2[ih][ia][ib]
                            *le2[ih][id][ic]
                            *le2[is][id][ic]
                        *ccc[ia][ic];
        }
        amp[is][ih] = dtemp;
    }
    return 0;

}




////////////////
int constructor_f_exciton_gij::GetFullRedAmplitudes11(double**** amp, double** cstrength)
{
    int num;
    evecs->GetSize(num,num);
    double** lev = evecs->data2D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        for( int id=0; id<num; id++ )
        {
            double dtemp=0;
            for( int sn=0; sn<num; sn++ )
            for( int sm=0; sm<num; sm++ )
            {
                dtemp += lev[ia][sn]
                            *lev[ib][sn]
                            *lev[ic][sm]
                            *lev[id][sm]
                        *cstrength[sn][sm];
            }
            amp[ia][ib][ic][id] = dtemp;
        }
    return 0;
}
int constructor_f_exciton_gij::GetFullRedAmplitudes22(double**** amp, double** cstrength)
{
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;
    
    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int ih=0; ih<nu2; ih++ )
    for( int id=0; id<nu2; id++ )
    for( int ij=0; ij<nu2; ij++ )
    { 
 
        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        for( int id=0; id<num; id++ )
        {
               dtemp +=      le2[is][ia][ib]
                            *le2[ih][ia][ib]
                            *le2[id][id][ic]
                            *le2[ij][id][ic]
                        *cstrength[ia][ic];
        }
        amp[is][ih][id][ij] = dtemp;
    }
    return 0;
}
int constructor_f_exciton_gij::GetFullRedAmplitudes21(double**** amp, double** cstrength)
{
    double** lev = evecs->data2D;
    int num;
    int nu2;
    evec2->GetSize(nu2,num,num);
    double*** le2 = evec2->data3D;

    // here we assume that all sites fluctuate with the same type of fluctuation independently.
    for( int is=0; is<nu2; is++ )
    for( int iz=0; iz<nu2; iz++ )
    for( int ie=0; ie<num; ie++ )
    for( int ix=0; ix<num; ix++ )
    {


        double dtemp=0;
        for( int ia=0; ia<num; ia++ )
        for( int ib=0; ib<num; ib++ )
        for( int ic=0; ic<num; ic++ )
        {
               dtemp +=      lev[ie][ia]
                            *lev[ix][ia]
                            *le2[is][ic][ib]
                            *le2[iz][ic][ib]
                        *cstrength[ia][ib];
        }

        amp[is][iz][ie][ix] = dtemp;
    }
    return 0;
}



int constructor_f_exciton_gij::GetFullAmplitudesMultimode(storage<double>& results, storage<double>& source)
{
    int howmuch = 1;
    int num,numq,numo;
    int nu2;
    double*** le2;
    if(evec2 != 0)
    {
        howmuch = 2;
        evec2->GetSize(nu2,num,num);
        le2 = evec2->data3D;
    }else{
        evecs->GetSize(num,num);
    }
    double** lev = evecs->data2D;
    double** cst = source.data2D;
    double*** cst2 = source.data3D;
    double***** amp = results.data5D;

    results.GetSize(numq,numq,numq,numq,numo);

	if(howmuch == 2)
	{
		if(numq != nu2+num+1)
    		{
        		std::cout<<"ERROR: the array sizes in constructor_f_exciton_gij::GetFullAmplitudes do not match\n";
		}
	}
	else if(numq != num+1)
    	{
        	std::cout<<"ERROR: the array sizes in constructor_f_exciton_gij::GetFullAmplitudes do not match\n";
    	}

	for( int io=0; io<numo; io++ ){

    if(source.CheckDimension()==2)
    {
        // uncorrelated
        // only 1 exc manifold
	{
	    // here we assume that all sites fluctuate with the same type of fluctuation independently.
	    for( int ia=0; ia<num; ia++ )
	        for( int ib=0; ib<num; ib++ )
	    for( int ic=0; ic<num; ic++ )
	        for( int id=0; id<num; id++ )
	        {
	            double dtemp=0;
	            for( int sn=0; sn<num; sn++ )
	            {
	                dtemp += lev[ia][sn]*lev[ib][sn]*lev[ic][sn]*lev[id][sn]*cst[io][sn];
	            }
	            amp[ia+1][ib+1][ic+1][id+1][io] += dtemp;
                }
        }
        if(howmuch == 2)  //  1 and 2 exc manifold
        {

		// 21
		    for( int is=0; is<nu2; is++ )
		    for( int iz=0; iz<nu2; iz++ )
		{
		    for( int ie=0; ie<num; ie++ )
		    for( int ix=0; ix<num; ix++ )
		    {


		        double dtemp=0;
		        for( int ia=0; ia<num; ia++ )
		        for( int ic=0; ic<num; ic++ )
		        {
		               dtemp +=  lev[ie][ia]
		                          *lev[ix][ia]
		                          *le2[is][ic][ia]
		                          *le2[iz][ic][ia]
		                        *cst[io][ia];
			}
			amp[is+1+num][iz+1+num][ie+1][ix+1][io] += dtemp;
			amp[ie+1][ix+1][is+1+num][iz+1+num][io] += dtemp;
		    }
                //22
   		    for( int ii=0; ii<nu2; ii++ )
    		    for( int ij=0; ij<nu2; ij++ )
    			{ 
 
		        double dtemp=0;
        		for( int ia=0; ia<num; ia++ )
        		for( int ib=0; ib<num; ib++ )
        		for( int id=0; id<num; id++ )
        		{
               dtemp += le2[is][ia][ib]
                            *le2[iz][ia][ib]
                            *le2[ii][id][ia]
                            *le2[ij][id][ia]*cst[io][ia];
        		}
			amp[is+1+num][iz+1+num][ii+1+num][ij+1+num][io] += dtemp;
			}
		}
	}
    }else
    {
        // correlated
        // only 1 exc manifold
	{
	    // here we assume that all sites fluctuate with the same type of fluctuation independently.
	    for( int ia=0; ia<num; ia++ )
	        for( int ib=0; ib<num; ib++ )
	    for( int ic=0; ic<num; ic++ )
	        for( int id=0; id<num; id++ )
	        {
	            double dtemp=0;
	            for( int sn=0; sn<num; sn++ )
	            for( int sm=0; sm<num; sm++ )
	            {
	                dtemp += lev[ia][sn]*lev[ib][sn]*lev[ic][sm]*lev[id][sm]*cst2[io][sn][sm];
	            }
	            amp[ia+1][ib+1][ic+1][id+1][io] += dtemp;
                }
        }
        if(howmuch == 2)  //  1 and 2 exc manifold
        {

		// 21
		    for( int is=0; is<nu2; is++ )
		    for( int iz=0; iz<nu2; iz++ )
		{
		    for( int ie=0; ie<num; ie++ )
		    for( int ix=0; ix<num; ix++ )
		    {


		        double dtemp=0;
		        for( int ia=0; ia<num; ia++ )
		        for( int ib=0; ib<num; ib++ )
		        for( int ic=0; ic<num; ic++ )
		        {
		               dtemp +=  lev[ie][ia]
		                          *lev[ix][ia]
		                          *le2[is][ic][ib]
		                          *le2[iz][ic][ib]
		                        *cst2[io][ia][ib];
			}
			amp[is+1+num][iz+1+num][ie+1][ix+1][io] += dtemp;
			amp[ie+1][ix+1][is+1+num][iz+1+num][io] += dtemp;
		    }
                //22
   		    for( int ii=0; ii<nu2; ii++ )
    		    for( int ij=0; ij<nu2; ij++ )
    			{ 
 
		        double dtemp=0;
        		for( int ia=0; ia<num; ia++ )
        		for( int ib=0; ib<num; ib++ )
        		for( int ic=0; ic<num; ic++ )
        		for( int id=0; id<num; id++ )
        		{
               dtemp += le2[is][ia][ib]
                            *le2[iz][ia][ib]
                            *le2[ii][id][ic]
                            *le2[ij][id][ic]*cst2[io][ia][ic];
        		}
			amp[is+1+num][iz+1+num][ii+1+num][ij+1+num][io] += dtemp;
			}
		}
	}

	}

		}

	return 0;
}

