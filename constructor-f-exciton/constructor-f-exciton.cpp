// everything is performed
#include"constructor-f-exciton.hpp"


#include"../eigen/eigen.hpp"
#include"../constants/constants.hpp"


constructor_f_exciton::constructor_f_exciton(
                      storage<double>& ham,
                      storage<dvector3d>& dips
                      )
{
    Init();
    
    int numl,numr;
    ham.GetSize(numl,numr);

    if(numl!=numr || numl==0 || numr == 0 )
    {
        cout<<"Error: wrong size of the Hamiltonian matrix\n";
    }
    
    hamloc = ham;
    dipsloc = dips;
}

constructor_f_exciton::constructor_f_exciton(
                                             int lev,
                                             storage<double>& ham,
                                             storage<dvector3d>& dips
                                             )
{

    if(!(lev == 2 || lev == 3))
    {
        cout<<"Error: The program supports only two- and three-level molecules\n";
        cout<<"Proceeding with three-level molecules\n";
    }
    Init();
    levels = lev;
    
    int numl,numr;
    ham.GetSize(numl,numr);

    if(numl!=numr || numl==0 || numr == 0 )
    {
        cout<<"Error: wrong size of the Hamiltonian matrix\n";
    }
    
    hamloc = ham;
    dipsloc = dips;
}


constructor_f_exciton::constructor_f_exciton(storage<double>& ham)
{
    cout<<"Warning: transition amplitudes shall not be used\n";
    
    Init();
    int numl,numr;
    
    ham.GetSize(numl,numr);

    if(numl!=numr || numl==0 || numr == 0 )
    {
        cout<<"Error: wrong size of the Hamiltonian matrix\n";
    }
    
    hamloc = ham;
}

constructor_f_exciton::constructor_f_exciton(int lev, storage<double>& ham)
{
    cout<<"Warning: transition amplitudes shall not be used\n";
    if(!(lev == 2 || lev == 3))
    {
        cout<<"Error: The program supports only two- and three-level molecules\n";
        cout<<"Warning: proceeding with three-level molecules\n";
    }
    Init();
    
    levels = lev;
    int numl,numr;
    
    ham.GetSize(numl,numr);

    if(numl!=numr || numl==0 || numr == 0 )
    {
        cout<<"Error: wrong size of the Hamiltonian matrix\n";
    }
    
    hamloc = ham;
}

constructor_f_exciton::constructor_f_exciton()
{
    cout<<"Error: one should not use \"empty\" constructor\n";
    Init();
}
constructor_f_exciton::~constructor_f_exciton()
{
}

int constructor_f_exciton::Init()
{
    // default value
    levels = 3;
        
    evals.SetDimension(1);
    evecs.SetDimension(2);
    edips.SetDimension(1);
    evals2.SetDimension(1);
    evecs2.SetDimension(2);
    edips2.SetDimension(2);
    
    hamloc.SetDimension(2);
    dipsloc.SetDimension(1);
    hamloc2.SetDimension(2);
    dipsloc2.SetDimension(2);
    
    hamK.SetDimension(2);
    dips2.SetDimension(1);

    lpos.SetDimension(1); // local coordinates of transitions
    lten.SetDimension(1); // local coordinates of transitions
    lmag.SetDimension(1); // local magnetic dipoles
    eten.SetDimension(1); // eigentensors of quadrupole transitions
    emag.SetDimension(1); // eigente magnetic transition dipoles
    lten2.SetDimension(2); // local coordinates of transitions
    lmag2.SetDimension(2); // local magnetic dipoles
    eten2.SetDimension(2); // eigentensors of quadrupole transitions
    emag2.SetDimension(2); // eigente magnetic transition dipoles

	manifold2sorter_u.SetDimension(1);
	manifold2sorter_l.SetDimension(1);

}

int constructor_f_exciton::AddKCouplings(storage<double>& iK)
{    
    int numl,numr;
    
    if(!iK.IsAlloc())
        return 2;
    
    if(iK.CheckDimension() == 1 && levels == 3)
    {
        
        // OK to proceed
        int numl;
        iK.GetSize(numl);// the size must be the same as of hamloc
        hamK.Allocate(numl,numl);
        for(int il=0; il<numl; il++)
            for(int ir=0; ir<numl; ir++)
            {
                if(il == ir)
                    hamK.data2D[il][ir]=iK.data1D[il];
                else
                    hamK.data2D[il][ir]=0.0;
            }
    }
    else if(iK.CheckDimension() == 2)
    {
        // OK to proceed
        int numl;
        iK.GetSize(numl,numl);// the size must be the same as of hamloc
        hamK.Allocate(numl,numl);
        for(int il=0; il<numl; il++)
            for(int ir=0; ir<numl; ir++)
            {
                    hamK.data2D[il][ir]=iK.data2D[il][ir];
            }
    }
    else
    {
        cout<<"Error: Dimensionality of K matrix is incompatible with the number of levels\n";
        return 1;
    }
}


int constructor_f_exciton::AddADipoles(storage<dvector3d>& id)
{    
    int numl,numr;
    
    if(levels == 3)
    {
        dips2=id;
    }
    else
    {
        cout<<"Error: anharmonic dipoles are incompatible with the number of levels\n";
        return 1;
    }
}




int constructor_f_exciton::MakeEigens1()
{
    if(evals.IsAlloc())
        return 0;
    
    int numl,numr;
    hamloc.GetSize(numl,numr);
    evals.Allocate(numl);
    evecs.Allocate(numl,numl);
    
    eigen obj(1);
    obj.GetEigenSolution(&hamloc, &evals, &evecs);
    return 0;
}

int constructor_f_exciton::MakeEigens2()
{
    if(evals2.IsAlloc())
        return 0;



    if(hamloc2.IsAlloc())
    {
 		int numl,numr;
 		hamloc2.GetSize(numl,numr);
 		evals2.Allocate(numl);
 		evecs2.Allocate(numl,numl);
     
 		eigen obj(1);
 		obj.GetEigenSolution(&hamloc2, &evals2, &evecs2);
 		return 0;
     }
     
     MakeHamLoc2();
     MakeEigens2();
}

int constructor_f_exciton::GetEvals(storage<double>& result)
{
    if(evals.IsAlloc())
    {
        result = evals;
        return 0;
    }
    
    MakeEigens1();
    GetEvals(result);
    
}

int constructor_f_exciton::GetEvecs(storage<double>& result)
{
    if(evecs.IsAlloc())
    {
        result = evecs;
        return 0;
    }
    
    MakeEigens1();
    GetEvecs(result);
    
}

int constructor_f_exciton::GetEvecsES(storage<double>& result)
{
    if(evecs.IsAlloc())
    {
        int numl;
        evals.GetSize(numl);
        for(int ie=0; ie<numl; ie++)
        for(int il=0; il<numl; il++)
            result.data2D[ie][il]=evecs.data2D[il][ie];
        return 0;
    }
    
    MakeEigens1();
    GetEvecsES(result);    
}

int constructor_f_exciton::GetEvals2(storage<double>& result)
{
    if(evals2.IsAlloc())
    {
        result = evals2;
        return 0;
    }
    
    MakeEigens2();
    GetEvals2(result);
    
}
void constructor_f_exciton::reformatF12()
{
    if(!evecs2.IsAlloc())
	return;

    if(evecs2.CheckDimension() ==2 )
    {
    // here I make the proper eigenvector representation
    // for 2-excitons
    
    // either two-level or three-level sites
    storage<double> et(3);
    
    int num1;
    evals.GetSize(num1);

    int sh = 0;
    
    int num2;
    evals2.GetSize(num2);

    if(levels == 2)
        // 2level
        sh = 0;
    else
        // 3level
        sh = num1;
    
    
    et.Allocate(num2,num1,num1);
    
    for(int ix = 0; ix<num2; ix++)
    {
	// overtones (for three level molecules)
	for(int ii=0;ii<sh;ii++)
        {
            et.data3D[ix][ii][ii] = evecs2.data2D[ii][ix];
        }
    
 	//Couplings
        for(int ii=0; ii<num1; ii++)
        for(int ij=0; ij<ii; ij++)
        {
            et.data3D[ix][ii][ij] = evecs2.data2D[(ii*ii-ii)/2+ij+sh][ix];
            et.data3D[ix][ij][ii] = et.data3D[ix][ii][ij];
        }
    }
    
    evecs2 = et;
    }
    
}
void constructor_f_exciton::reformatF21()
{
    if(!evecs2.IsAlloc())
	return;

    if(evecs2.CheckDimension() ==3 )
    {
    
    // either two-level or three-level sites
    storage<double> et(2);
    
    int num1;
    evals.GetSize(num1);

    int sh = 0;
    
    
    int num2;
    evals2.GetSize(num2);

    if(levels == 2)
        // 2level
        sh = 0;
    else
        // 3level
        sh = num1;
    
    
    et.Allocate(num2,num2);
    
    for(int ix = 0; ix<num2; ix++)
    {
	// overtones (for three level molecules)
	for(int ii=0;ii<sh;ii++)
        {
            et.data2D[ii][ix] = evecs2.data3D[ix][ii][ii];
        }
    
 	//Couplings
        for(int ii=0; ii<num1; ii++)
        for(int ij=0; ij<ii; ij++)
        {
            et.data2D[(ii*ii-ii)/2+ij+sh][ix] = evecs2.data3D[ix][ii][ij];
        }
    }
    
    evecs2 = et;
    }

}
int constructor_f_exciton::GetEvecs2(storage<double>& result)
{
    if(evecs2.IsAlloc() != 0)
    {
        result = evecs2;
        return 0;
    }
    
    MakeEigens2();
    GetEvecs2(result);
    
}
int constructor_f_exciton::GetStateEvecs2()
{
	if(!evecs2.IsAlloc())
		return 0;
	else if(evecs2.CheckDimension()==2)
		return 1;

	else return 2;
}


int constructor_f_exciton::GetHam2(storage<double>& result)
{
    if(hamloc2.IsAlloc())
    {
        result = hamloc2;
        return 0;
    }
    MakeHamLoc2();
    GetHam2(result);    
}

int constructor_f_exciton::GetHamiltonian2LocalSorters(storage<int>& sl,storage<int>&sr)
{
	if(!hamloc2.IsAlloc())
	{
		MakeHamLoc2();
	}
	sl = manifold2sorter_u;
	sr = manifold2sorter_l;
	return 0;
}

int constructor_f_exciton::MakeHamLoc2()
{
    // this constructs the two-exciton Hamiltonian for 
    // either two-level or three-level sites
    if(hamloc2.IsAlloc())
        return 0;
    
    int num1;
    hamloc.GetSize(num1,num1);
    int sh = 0;
    
    
    int num2 = (num1*num1-num1)/2;

    if(levels==3)//three-level systems
	{
		sh = num1;
        num2 = (num1*num1+num1)/2;
	}
    


    hamloc2.Allocate(num2,num2);
    manifold2sorter_u.Allocate(num2);
    manifold2sorter_l.Allocate(num2);
    
    double** lhamloc = hamloc.data2D;
    double** lhamloc2 = hamloc2.data2D;
    
    double** lhamK = 0;
    if(hamK.IsAlloc())
        lhamK = hamK.data2D;

	int absolutestate=0;
    
	// overtones (for three level molecules)
	for(int ii=0;ii<sh;ii++)
    {
        if(hamK.IsAlloc())
            lhamloc2[ii][ii] = 2*lhamloc[ii][ii] + lhamK[ii][ii];
        else
            lhamloc2[ii][ii] = 2*lhamloc[ii][ii];

	manifold2sorter_u.data1D[absolutestate]=ii;
	manifold2sorter_l.data1D[absolutestate]=ii;
	absolutestate++;
    }
    
    // Energies of two-excitons
    for(int ii=0;ii<num1; ii++)
        for(int jj=0;jj<ii;jj++)
        {
            int kk =(ii*ii-ii)/2+jj+sh;
            
            if(hamK.IsAlloc())
                lhamloc2[kk][kk] = lhamloc[ii][ii]+lhamloc[jj][jj]
                    + lhamK[ii][jj];
            else
                lhamloc2[kk][kk] = lhamloc[ii][ii]+lhamloc[jj][jj];

	manifold2sorter_u.data1D[absolutestate]=ii;
	manifold2sorter_l.data1D[absolutestate]=jj;
	absolutestate++;
        }
    
    
 	//Couplings
    for(int ii=0; ii<num1; ii++)
        for(int ij=0; ij<ii; ij++)
        {
            
        
            for(int ik=0; ik<num1; ik++)
                for(int il=0; il<ik; il++)
                {
                    if((ii==ik)&&(ij!=il))
                    {
                        lhamloc2[(ii*ii-ii)/2+ij+sh][(ik*ik-ik)/2+il+sh] = lhamloc[ij][il];
                        lhamloc2[(ik*ik-ik)/2+il+sh][(ii*ii-ii)/2+ij+sh] = lhamloc[ij][il];
                    }
                    
                    if((ij==ik)&&(ii!=il))
                    {
                        lhamloc2[(ii*ii-ii)/2+ij+sh][(ik*ik-ik)/2+il+sh] = lhamloc[ii][il];
                        lhamloc2[(ik*ik-ik)/2+il+sh][(ii*ii-ii)/2+ij+sh] = lhamloc[ii][il];;
                    }
                    
                    if((ii==il)&&(ij!=ik))
                    {
                        lhamloc2[(ii*ii-ii)/2+ij+sh][(ik*ik-ik)/2+il+sh] = lhamloc[ij][ik];
                        lhamloc2[(ik*ik-ik)/2+il+sh][(ii*ii-ii)/2+ij+sh] = lhamloc[ij][ik];
                    }
                    
                    if((ij==il)&&(ii!=ik))
                    {
                        lhamloc2[(ii*ii-ii)/2+ij+sh][(ik*ik-ik)/2+il+sh] = lhamloc[ii][ik];
                        lhamloc2[(ik*ik-ik)/2+il+sh][(ii*ii-ii)/2+ij+sh] = lhamloc[ii][ik];
                    }
                    
                }
        }
            
    // combination band and overtone couplings (only for three-level systems)
    bool dorescaling = false;
    if(dipsloc.IsAlloc() && dips2.IsAlloc())
    {
        cout<<"# constructor_f_exciton::MakeHamLoc2(): rescaling overtone-combination couplings based on dipole strengths\n";
        dorescaling = true;
    }
    
    for(int ii=0;   ii<sh; ii++)
        for(int ij=0;   ij<ii; ij++)
        {
            
            double scale_n = constants::root2;
            double scale_d = 1.0;
            
            if(dorescaling)
            {
                scale_d = dipsloc.data1D[ii].Amplitude()*dipsloc.data1D[ij].Amplitude();
            }
            
            //double coupl = constants::root2*lhamloc[ii][ij];
                    
            // ii overtone
            if(dorescaling)
            {
                scale_n = constants::root2
                        *dvector3d(dipsloc.data1D[ii]+dips2.data1D[ii]).Amplitude()
                        *dipsloc.data1D[ij].Amplitude();
            }
            
            lhamloc2[ii][(ii*ii-ii)/2+ij+sh] = scale_n/scale_d*lhamloc[ii][ij];
            lhamloc2[(ii*ii-ii)/2+ij+sh][ii] = scale_n/scale_d*lhamloc[ii][ij];
                    
            // ij overtone
            if(dorescaling)
            {
                scale_n = constants::root2
                        *dvector3d(dipsloc.data1D[ij]+dips2.data1D[ij]).Amplitude()
                        *dipsloc.data1D[ii].Amplitude();
            }
            lhamloc2[ij][(ii*ii-ii)/2+ij+sh] = scale_n/scale_d*lhamloc[ii][ij];
            lhamloc2[(ii*ii-ii)/2+ij+sh][ij] = scale_n/scale_d*lhamloc[ii][ij];
        }
            
    
    //  that is it
}

int constructor_f_exciton::MakeEdips1()
{
    if(edips.IsAlloc())
        return 0;
    
    if(!dipsloc.IsAlloc())
    {
        cout<<"Error: Local dipoles are absent: cannot make eigendipoles!\n";
        return 1;
        
    }
    
    if(evals.IsAlloc())
    {
        int numl;
        evals.GetSize(numl);
        
        edips.Allocate(numl);
        for(int ind = 0; ind< numl; ind ++)
        {
            edips.data1D[ind]=0.;
            for(int is = 0; is< numl; is++)
                edips.data1D[ind] += dipsloc.data1D[is]*evecs.data2D[is][ind];
        }

	if(lpos.IsAlloc()) // making transition tensors
	{
        eten.Allocate(numl);
        for(int ind = 0; ind< numl; ind ++)
        {
            eten.data1D[ind]=0.;
            for(int is = 0; is< numl; is++)
                eten.data1D[ind] += dtensor3x3(lpos.data1D[is],dipsloc.data1D[is])*evecs.data2D[is][ind];
        }
	}

	if(lmag.IsAlloc())
	{
        emag.Allocate(numl);
        for(int ind = 0; ind< numl; ind ++)
        {
            emag.data1D[ind]=0.;
            for(int is = 0; is< numl; is++)
                emag.data1D[ind] += lmag.data1D[is]*evecs.data2D[is][ind];
        }
	}

        return 0;
    }
    
    MakeEigens1();
    MakeEdips1();
    return 0;
}
int constructor_f_exciton::Makedipsloc2()
{
    // this constructs the two-exciton Hamiltonian for 
    // either two-level or three-level sites
    if(dipsloc2.IsAlloc())
        return 0;
    
    if(!dipsloc.IsAlloc())
    {
        cout<<"Error: Local dipoles are absent: cannot make eigendipoles!\n";
        return 1;
        
    }
    
    int num1;
    hamloc.GetSize(num1,num1);
    int sh = 0;
    int num2 = (num1*num1-num1)/2;
    
    if(levels==3)//three-level systems
	{
	sh = num1;
        num2 = (num1*num1+num1)/2;
	}
    
    dipsloc2.Allocate(num1,num2);
	if(lpos.IsAlloc())
		lten2.Allocate(num1,num2);
	if(lmag.IsAlloc())
		lmag2.Allocate(num1,num2);

    constants cst;
    
    // overtones (for three level molecules)
    if(dips2.IsAlloc()){
    for(int ii=0;ii<sh;ii++)
    {
        dipsloc2.data2D[ii][ii] = dipsloc.data1D[ii]+dips2.data1D[ii];
        dipsloc2.data2D[ii][ii] = dipsloc2.data2D[ii][ii]*cst.root2;
    }}
    else{
    for(int ii=0;ii<sh;ii++)
    {
        dipsloc2.data2D[ii][ii] = dipsloc.data1D[ii]*cst.root2;
    }}
        
    
    // dipoles  of two-excitons combination bands
    for(int ii=0;ii<num1; ii++)
        for(int jj=0;jj<ii;jj++)
        {
            int kk =(ii*ii-ii)/2+jj+sh;            
            dipsloc2.data2D[ii][kk]= dipsloc.data1D[jj];
            dipsloc2.data2D[jj][kk]= dipsloc.data1D[ii];
        }
    
	if(lpos.IsAlloc())
	{
	    // overtones (for three level molecules)
	    for(int ii=0;ii<sh;ii++)
	    {
	        lten2.data2D[ii][ii] = dtensor3x3(lpos.data1D[ii],dipsloc.data1D[ii]*cst.root2);
	    }
        
    	    // dipoles  of two-excitons combination bands
		for(int ii=0;ii<num1; ii++)
    		    for(int jj=0;jj<ii;jj++)
		{
		int kk =(ii*ii-ii)/2+jj+sh;            
		lten2.data2D[ii][kk]= dtensor3x3(lpos.data1D[jj],dipsloc.data1D[jj]);
		lten2.data2D[jj][kk]= dtensor3x3(lpos.data1D[ii],dipsloc.data1D[ii]);
		}
	}
	if(lmag.IsAlloc())
	{
	    // overtones (for three level molecules)
	    for(int ii=0;ii<sh;ii++)
	    {
	        lmag2.data2D[ii][ii] = lmag.data1D[ii]*cst.root2;
	    }
        
    	    // dipoles  of two-excitons combination bands
		for(int ii=0;ii<num1; ii++)
    		    for(int jj=0;jj<ii;jj++)
		{
		int kk =(ii*ii-ii)/2+jj+sh;            
		lmag2.data2D[ii][kk]= lmag.data1D[jj];
		lmag2.data2D[jj][kk]= lmag.data1D[ii];
		}
	}
    //  that is it
}

int constructor_f_exciton::MakeEdips2()
{
    if(edips2.IsAlloc())
        return 0;
    
    if(dipsloc.IsAlloc() == 0 )
    {
        cout<<"Error: Local dipoles are absent: cannot make eigendipoles!\n";
        return 1;
    }
    
    if(!dipsloc2.IsAlloc())
    {
        Makedipsloc2();
    }

    if(!evals.IsAlloc())
    {
        MakeEigens1();
    }

    if(!evals2.IsAlloc())
    {
        MakeEigens2();
    }
    
    int saveevecs2ready = evecs2.CheckDimension();
    if(saveevecs2ready == 3)
        reformatF21();

    
    int num1;
    evals.GetSize(num1);

    int num2;
    evals2.GetSize(num2);
        
    edips2.Allocate(num1,num2); // (one - to - two) transition
    for(int i1 = 0; i1< num1; i1 ++)
    for(int i2 = 0; i2< num2; i2 ++)
    {
        edips2.data2D[i1][i2]=0.;

        for(int is1 = 0; is1< num1; is1++)
        for(int is2 = 0; is2< num2; is2++)
            edips2.data2D[i1][i2] += (dipsloc2.data2D[is1][is2]
                    *(evecs.data2D[is1][i1]
                      *evecs2.data2D[is2][i2]));
    }

	if(lpos.IsAlloc())
	{
	    eten2.Allocate(num1,num2); // (one - to - two) transition
	    for(int i1 = 0; i1< num1; i1 ++)
	    for(int i2 = 0; i2< num2; i2 ++)
	    {
	        eten2.data2D[i1][i2]=0.;
	
	        for(int is1 = 0; is1< num1; is1++)
	        for(int is2 = 0; is2< num2; is2++)
	            eten2.data2D[i1][i2] += (lten2.data2D[is1][is2]
	                    *(evecs.data2D[is1][i1]
	                      *evecs2.data2D[is2][i2]));
	    }
	}
	if(lmag.IsAlloc())
	{
	    emag2.Allocate(num1,num2); // (one - to - two) transition
	    for(int i1 = 0; i1< num1; i1 ++)
	    for(int i2 = 0; i2< num2; i2 ++)
	    {
	        emag2.data2D[i1][i2]=0.;
	
	        for(int is1 = 0; is1< num1; is1++)
	        for(int is2 = 0; is2< num2; is2++)
	            emag2.data2D[i1][i2] += (lmag2.data2D[is1][is2]
	                    *(evecs.data2D[is1][i1]
	                      *evecs2.data2D[is2][i2]));
	    }
	}
    
    if(saveevecs2ready == 3)
        reformatF12();

    return 0;
}

int constructor_f_exciton::AddDipolePositions(storage<dvector3d>& ip)
{
	lpos = ip;
	return 0;
}
int constructor_f_exciton::AddMagneticDipoles(storage<dvector3d>& im)
{
	lmag = im;
	return 0;
}

int constructor_f_exciton::GetETens(storage<dtensor3x3>& ten)
{
    if(edips.IsAlloc())
    {
	ten = eten;
	return 0;
    }
    
    if(MakeEdips1() == 0);
        return GetETens(ten);
}
int constructor_f_exciton::GetEMags(storage<dvector3d>& mag)
{
    if(edips.IsAlloc())
    {
	mag = emag;
	return 0;
    }
    
    if(MakeEdips1() == 0);
        return GetEMags(mag);
}
int constructor_f_exciton::GetETens2(storage<dtensor3x3>& ten)
{
    if(edips2.IsAlloc())
    {
	ten = eten2;
	return 0;
    }
    
    MakeEdips2();
    return GetETens2(ten);    
}
int constructor_f_exciton::GetEMags2(storage<dvector3d>& mag)
{
    if(edips2.IsAlloc())
    {
	mag = emag2;
	return 0;
    }
    
    MakeEdips2();
    return GetEMags2(mag);    
}




int constructor_f_exciton::GetEdips(storage<dvector3d>& result)
{
    if(edips.IsAlloc())
    {
        result = edips;
        return 0;
    }
    
    if(MakeEdips1() == 0);
        return GetEdips(result);

}
int constructor_f_exciton::GetEdips2(storage<dvector3d>& result)
{
    if(edips2.IsAlloc())
    {
        result = edips2;
        return 0;
    }
    
    MakeEdips2();
    return GetEdips2(result);    
}

int constructor_f_exciton::GetDips2(storage<dvector3d>& result)
{
    if(dipsloc2.IsAlloc())
    {
        result = dipsloc2;
        return 0;
    }
    
    Makedipsloc2();
    return GetDips2(result);
}
int constructor_f_exciton::GetDips(storage<dvector3d>& result)
{
    if(dipsloc.IsAlloc())
    {
        result = dipsloc;
        return 0;
    }
	std::cout<<"Error: local dipoles are missing in constructor_f_exciton::GetDips()\n";
	return 1;
}

int constructor_f_exciton::GetTens(storage<dtensor3x3>& ten)
{
	if(!lten.IsAlloc())
		Makedipsloc();
	ten = lten;

	if(lten.IsAlloc())
		return 0;
	else return 1;
}
int constructor_f_exciton::GetMags(storage<dvector3d>& mag)
{
	mag= lmag;

	if(lmag.IsAlloc())
		return 0;
	else return 1;
}
int constructor_f_exciton::GetTens2(storage<dtensor3x3>& ten)
{
	if(!dipsloc2.IsAlloc())
		Makedipsloc2();
	ten= lten2;

	if(lten2.IsAlloc())
		return 0;
	else return 1;
}
int constructor_f_exciton::GetMags2(storage<dvector3d>& mag)
{
	if(!dipsloc2.IsAlloc())
		Makedipsloc2();
	mag= lmag2;

	if(lmag2.IsAlloc())
		return 0;
	else return 1;
}



int constructor_f_exciton::Makedipsloc()
{
    if(dipsloc.IsAlloc())
    {
        // local dipoles are present. OK have to make tensors 
	if(lpos.IsAlloc()) // making transition tensors
	{
        lten.Allocate(lpos.GetSize());
        for(int ind = 0; ind< lpos.GetSize(); ind ++)
        {
            lten.data1D[ind]= dtensor3x3(lpos.data1D[ind],dipsloc.data1D[ind]);
        }
	}

        return 0;
    }
	std::cout<<"Error: local dipoles are missing in constructor_f_exciton::Makedipsloc()\n";
	return 1;

}
