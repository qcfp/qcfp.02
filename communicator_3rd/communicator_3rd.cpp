#include"../toolsFFT/toolsFFT.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../storage/storage.hpp"
#include"../constants/constants.hpp"
#include"../feinman2sideddiagram1/feinman2sideddiagram1.hpp"
#include"../interpolationF2d/interpolationF2d.hpp"
#include"../feinman2sideddiagram3/feinman2sideddiagram3.hpp"
#include"../propagatorM/propagatorM.hpp"
#include"../numericalSD/numericalSD.hpp"
#include"../calculator_redfield/calculator_redfield.hpp"
#include"../constructor-f-exciton/constructor-f-exciton.hpp"
#include"../constructor-f-exciton-gij/constructor-f-exciton-gij.hpp"

#include"communicator_3rd.hpp"

communicator_3rd::communicator_3rd()
{
    Setup();

}

void communicator_3rd::MakeESSystem()
{
	MakeESSystem(0);
}

void communicator_3rd::MakeESSystem(int flag)
{
    
    int numBosc = mfun.GetSize();
    
    // going to eigensystem
    // system in eigenstate basis
    storage<double> leval1(1);
    storage<double> levecs(2);
    storage<dvector3d> edip1(1);

    storage<double> leval2(1);
    storage<double> levec2(3); //[ix][ia][ib]
    storage<dvector3d> edip2(2);


    leval1.Allocate(numE);
    levecs.Allocate(numE,numE);
    edip1.Allocate(numE);

    if(band == 2)
    {
      leval2.Allocate(numF);
      levec2.Allocate(numF,numE,numE);
      edip2.Allocate(numE,numF);
    }

    // now making the two-level systems
    int numlev = 2;
    if(bosonic) numlev = 3;
    constructor_f_exciton exc_system(numlev,ham,dip);

    
    // reorganization energies
    cout<<"# Shifting electronic Hamiltonian energies up by site reorganization energies.\n";
    cout<<"# Site reorganization energies:\n";
    for(int indl=0; indl<numE; indl++)
    {
        double reorganization = 0.0;
        
        // 	if(flag>0) //  shared baths approach where on the input are directly coupling constants
        // 	{
        //                 // first I calculate all site reorganization energies
        //                 if(coupl.CheckDimension()==2)
        //                 {
        //                     for(int ind = 0; ind < numBosc; ind++)
        //                     reorganization += -pow(coupl.data2D[ind][indl],2)*(mfun.data1D[ind].Get(0.0)).imag();
        //                 }
        //                 else
        //                 {
        // 			std::cout<<"Error in communicator_3rd::MakeESSystem: this approach is not possible\n";
        //                 }
        //                 
        // 	    if((band == 2) && bosonic) 
        //     		{
        //                 	anh.data2D[indl][indl] += reorganization*2;
        //     		}
        //                 
        //                 cout<<reorganization<<"\n";
        //             }
        // 
        // 
        // 	}
        // 	else
        // 	{
        // traditional approach
        // first I calculate all site reorganization energies
        if(coupl.CheckDimension()==2)
        {
            for(int ind = 0; ind < numBosc; ind++)
                reorganization += -coupl.data2D[ind][indl]*(mfun.data1D[ind].Get(0.0)).imag();
        }
        else
        {
            for(int ind = 0; ind < numBosc; ind++)
                reorganization += -coupl.data3D[ind][indl][indl]*(mfun.data1D[ind].Get(0.0)).imag();
        }
        
        
        
        ham.data2D[indl][indl] += reorganization;
        if((band == 2) && bosonic) 
        {
            //            anh.data2D[indl][indl] += reorganization*2;
            anh.data2D[indl][indl] += reorganization*4;
        }
        cout<<indl<<" "<<reorganization<<"\n";
    }
    
    std::cout<<" Hamiltonian after reorganizations\n";
    for(int indl=0; indl<numE; indl++)
        for(int indr=0; indr<=indl; indr++)
        {
            cout<<ham.data2D[indl][indr];
            if(indl==indr)cout<<"\n"; else cout<<"\t";
        }
        
    if(band == 2 && anh.IsAlloc()) 
    std::cout<<" Anharmonicities after reorganizations\n";
    if(band == 2 && anh.IsAlloc()) 
    for(int indl=0; indl<numE; indl++)
    for(int indr=0; indr<=indl; indr++)
    {
        cout<<anh.data2D[indl][indr];
        if(indl==indr)cout<<"\n"; else cout<<"\t";
    }

    if(band == 2) exc_system.AddKCouplings(anh);
    if(band == 2 && bosonic && dip2Corrections.IsAlloc()) 
        exc_system.AddADipoles(dip2Corrections);

	if(pos.IsAlloc())
		exc_system.AddDipolePositions(pos);
	if(d_m.IsAlloc())
		exc_system.AddMagneticDipoles(d_m);

    exc_system.GetEvals(leval1);
    exc_system.GetEvecsES(levecs);
    exc_system.GetEdips(edip1);
    if(band == 2) exc_system.GetEvals2(leval2);
    if(band == 2) exc_system.reformatF12();
    if(band == 2) exc_system.GetEvecs2(levec2);
    if(band == 2) exc_system.GetEdips2(edip2);
    if(band == 2) exc_system.GetHamiltonian2LocalSorters(ham2sorter_l,ham2sorter_r);




	storage<dtensor3x3> leten1(1);
	storage<dtensor3x3> leten2(2);
	storage<dvector3d> lemag1(1);
	storage<dvector3d> lemag2(2);
	if(pos.IsAlloc())
	{
		exc_system.GetETens(leten1);
		if(band == 2) exc_system.GetETens2(leten2);
	}
	if(d_m.IsAlloc())
	{
		exc_system.GetEMags(lemag1);
		if(band == 2) exc_system.GetEMags2(lemag2);
	}

    cout<<"#########################################\n";
    cout<<"# Eigensystem description:\n";
    cout<<"# number of zero-band levels:\n";
    cout<<1<<"\n";
    cout<<"# number of one-band levels:\n";
    cout<<numE<<"\n";
    if(band == 2) cout<<"# number of two-band levels:\n";
    if(band == 2) cout<<numF<<"\n";

    // print eigenvalues
    cout<<"# zero-band level energies:\n";
    cout<<0.0<<"\n";

    // print eigenvalues
    cout<<"# one-band level energies:\n";
    for(int ind=0; ind<numE; ind++)
        cout<<leval1.data1D[ind]<<"\n";

    // print eigenvalues 2
    if(band == 2) cout<<"# two-band level energies:\n";
    if(band == 2) 
      for(int ind=0; ind<numF; ind++)
        cout<<leval2.data1D[ind]<<"\n";


    // print eigenvectors
    cout<<"# Evecs SE (eigenvectors in columns):\n";
    for(int indl=0; indl<numE; indl++)
    {
        cout<<"# ";
        for(int inde=0; inde<numE; inde++)
        {
            cout<<levecs.data2D[inde][indl];
            if(inde == numE-1) cout<<"\n";
            else  cout<<"\t";
        }
    }


    cout<<"# Participation (delocalization) numbers for E:\n";
    for(int inde=0; inde<numE; inde++)
    {
        double delocalization = 0;
        
        for(int indl=0; indl<numE; indl++)
        {
            delocalization += (levecs.data2D[inde][indl]*levecs.data2D[inde][indl]*levecs.data2D[inde][indl]*levecs.data2D[inde][indl]);
        }
        cout<<1.0/delocalization<<"\n";
    }


    // print eigenvectors F
    if(band == 2) cout<<"# Evecs SSF (eigenvectors in columns):\n";
    if(band == 2) 
    for(int indl=0; indl<numE; indl++)
    {
        cout<<"# ";
        for(int indr=0; indr<=indl; indr++)
        {
            cout<<indl<<"-"<<indr<<"\t";
            for(int inde=0; inde<numF; inde++)
            {
                cout<<levec2.data3D[inde][indl][indr];
                if(inde == numF-1) cout<<"\n";
                else  cout<<"\t";
            }

        }
    }

    if(band == 2) cout<<"# Participation (delocalization) numbers for F in real space:\n";
    if(band == 2) 
    for(int inde=0; inde<numF; inde++)
    {
        double delocalization = 0;
        
        for(int indl=0; indl<numE; indl++)
        for(int indj=0; indj<=indl; indj++)
        {
            delocalization += levec2.data3D[inde][indl][indj]*levec2.data3D[inde][indl][indj]*levec2.data3D[inde][indl][indj]*levec2.data3D[inde][indl][indj];
        }
        cout<<1.0/delocalization<<"\n";
    }

    
    cout<<"# Transition dipoles:\n";
    cout<<"# number of enries:\n";
    cout<<numE+numE*numF<<"\n";


    // print eigendipoles
    cout<<"# Edips 0-E:\n";
    for(int ind=0; ind<numE; ind++)
    {
        cout<<"0\t"<<ind+1<<"\t";
        cout<<edip1.data1D[ind].x()<<"\t";
        cout<<edip1.data1D[ind].y()<<"\t";
        cout<<edip1.data1D[ind].z()<<"\n";
    }

    // print eigendipoles
    if(band == 2) cout<<"# Edips E-F:\n";
    if(band == 2) 
    for(int ind=0; ind<numE; ind++)
        for(int inf=0; inf<numF; inf++)
        {
            cout<<ind+1<<"\t"<<inf+1+numE<<"\t";
            //        cout<<ind<<"\t"<<inf<<"\t";
            cout<<edip2.data2D[ind][inf].x()<<"\t";
            cout<<edip2.data2D[ind][inf].y()<<"\t";
            cout<<edip2.data2D[ind][inf].z()<<"\n";
        }

    /*
     cout<<"all transition dipoles:\n";
     cout<<numE+numE*numF<<"\n";
     for(int ind=0; ind<numE; ind++)
     {
     cout<<0<<"\t"<<ind+1<<"\t";
     cout<<edips.data1D[ind].x()<<"\t";
     cout<<edips.data1D[ind].y()<<"\t";
     cout<<edips.data1D[ind].z()<<"\n";
     }
     for(int ind=0; ind<numE; ind++)
     for(int inf=0; inf<numF; inf++)
     {
     cout<<ind+1<<"\t"<<inf+1+numE<<"\t";
     cout<<edip2.data2D[ind][inf].x()<<"\t";
     cout<<edip2.data2D[ind][inf].y()<<"\t";
     cout<<edip2.data2D[ind][inf].z()<<"\n";
     }

     */



    // translating to eigenstate system

    evec0.Allocate(1,1);
    evec0.data2D[0][0]=1;
    evec1 = levecs;
    if(band == 2) evec2 = levec2;

	if(flag>0) //  shared baths approach where on the input are directly coupling constants
	{
		// do not transform bath since the computation will take place in the site basis
	}
	else
	{
		// traditional approach

    // fluctuation amplitudes
    gij.Allocate(numT,numT,numBosc);
    kij.Allocate(numT,numT,numBosc);



    cout<<"# making spectral shapes:\n";

    storage<double>* levec2marker=0;
    if(band == 2) levec2marker = &levec2;

    constructor_f_exciton_gij exc_fluct(&levecs,levec2marker);

    for(int ind = 0; ind < numBosc; ind++)
    {

        // making diagonal and off-diagonal fluctuation amplitudes
        storage<double> gamp(2);

        // gij
        gamp.Allocate(numE,numE);

        if(coupl.CheckDimension()==2)
            exc_fluct.GetGijAmplitudes11(&gamp,coupl.data2D[ind]);
        else
            exc_fluct.GetGijAmplitudesC11(gamp,coupl.data3D[ind]);
            
        for(int ia=0;ia<numE;ia++)
            for(int ie=0;ie<numE;ie++)
                gij.data3D[ia+1][ie+1][ind] += gamp.data2D[ia][ie];
        gamp.Delete();

        if(band == 2) 
        {
        gamp.Allocate(numF,numE);

        if(coupl.CheckDimension()==2)
            exc_fluct.GetGijAmplitudes21(&gamp,coupl.data2D[ind]);
        else
            exc_fluct.GetGijAmplitudesC21(gamp,coupl.data3D[ind]);

        for(int ia=0;ia<numF;ia++)
            for(int ie=0;ie<numE;ie++)
            {
                gij.data3D[ia+1+numE][ie+1][ind] += gamp.data2D[ia][ie];
                gij.data3D[ie+1][ia+1+numE][ind] += gamp.data2D[ia][ie];
            }
        gamp.Delete();

        gamp.Allocate(numF,numF);
        
        if(coupl.CheckDimension()==2)
            exc_fluct.GetGijAmplitudes22(&gamp,coupl.data2D[ind]);
        else
            exc_fluct.GetGijAmplitudesC22(gamp,coupl.data3D[ind]);

        for(int ia=0;ia<numF;ia++)
            for(int ie=0;ie<numF;ie++)
                gij.data3D[ia+1+numE][ie+1+numE][ind] += gamp.data2D[ia][ie];
        gamp.Delete();
        }
        // kij
        gamp.Allocate(numE,numE);
        
        if(coupl.CheckDimension()==2)
            exc_fluct.GetRedAmplitudes11(&gamp,coupl.data2D[ind]);
        else
            exc_fluct.GetRedAmplitudesC11(gamp,coupl.data3D[ind]);

        for(int ia=0;ia<numE;ia++)
            for(int ie=0;ie<numE;ie++)
                kij.data3D[ia+1][ie+1][ind] += gamp.data2D[ia][ie];
        gamp.Delete();

        if(band == 2) 
        {
        gamp.Allocate(numF,numF);
        
        if(coupl.CheckDimension()==2)
            exc_fluct.GetRedAmplitudes22(&gamp,coupl.data2D[ind]);
        else
            exc_fluct.GetRedAmplitudesC22(gamp,coupl.data3D[ind]);

        for(int ia=0;ia<numF;ia++)
            for(int ie=0;ie<numF;ie++)
                kij.data3D[ia+1+numE][ie+1+numE][ind] += gamp.data2D[ia][ie];
        gamp.Delete();
        }

    }

//    if(nonsecular)
//    {
	    if(!zijkl.IsAlloc())
		    zijkl.Allocate(numG+numE+numF,numG+numE+numF,numG+numE+numF,numG+numE+numF,numBosc);
	    exc_fluct.GetFullAmplitudesMultimode(zijkl, coupl);
//    }

    for(int ind = 0; ind < numBosc; ind++)
    {
        cout<<"# offdiagonal s-b coupling correlation triangular matrix\n";
        for(int ia = 0; ia<numT; ia++)
            for(int ib = 0; ib<=ia; ib++)
            {
                cout<< kij.data3D[ia][ib][ind];

                if(ib == ia) cout<<"\n";
                else cout<<"\t";
            }
    

        cout<<"# diagonal  s-b coupling correlation triangular matrix\n";
        for(int ia = 0; ia<numT; ia++)
            for(int ib = 0; ib<=ia; ib++)
            {
                cout<< gij.data3D[ia][ib][ind];

                if(ib == ia) cout<<"\n";
                else cout<<"\t";
            }
    }


		}/// closing transformation of bath

    ////////////////////////////
    // eigenlevel total system:


    evals.Allocate(numT);
    edips.Allocate(numT,numT);

	if(pos.IsAlloc())
		eten.Allocate(numT,numT);
	if(d_m.IsAlloc())
		ed_m.Allocate(numT,numT);

    for(int ia = 0; ia<numE; ia++)
    {
        double reor = 0.0;
        for(int ib=0;ib<numBosc;ib++)
            reor += -(gij.data3D[ia+1][ia+1][ib]*mfun.data1D[ib].Get(0.0)).imag();
        // putting only electronic energies
        // these later will be shifted by reorganization energies again
        evals.data1D[ia+1]=leval1.data1D[ia]-reor;
        edips.data2D[0][ia+1]=edip1.data1D[ia];
        edips.data2D[ia+1][0]=edip1.data1D[ia];

	if(pos.IsAlloc())
	{
	        eten.data2D[0][ia+1]=leten1.data1D[ia];
	        eten.data2D[ia+1][0]=leten1.data1D[ia];
	}
	if(d_m.IsAlloc())
	{
	        ed_m.data2D[0][ia+1]=lemag1.data1D[ia];
	        ed_m.data2D[ia+1][0]=lemag1.data1D[ia];
	}
    }
    if(band == 2) 
    for(int ia = 0; ia<numF; ia++)
    {
        double reor = 0.0;
        for(int ib=0;ib<numBosc;ib++)
            reor += -(gij.data3D[ia+1+numE][ia+1+numE][ib]*mfun.data1D[ib].Get(0.0)).imag();
        // putting only electronic energies
        // these later will be shifted by reorganization energies again

        evals.data1D[ia+1+numE]=leval2.data1D[ia]-reor;
        for(int ib = 0; ib<numE; ib++)
        {
            edips.data2D[ia+1+numE][ib+1]=edip2.data2D[ib][ia];
            edips.data2D[ib+1][ia+1+numE]=edip2.data2D[ib][ia];


		if(pos.IsAlloc())
		{
	            eten.data2D[ia+1+numE][ib+1]=leten2.data2D[ib][ia];
	            eten.data2D[ib+1][ia+1+numE]=leten2.data2D[ib][ia];
		}
		if(d_m.IsAlloc())
		{
	            ed_m.data2D[ia+1+numE][ib+1]=lemag2.data2D[ib][ia];
	            ed_m.data2D[ib+1][ia+1+numE]=lemag2.data2D[ib][ia];
		}

        }
    }

    // population propagation time
    //cout<<"# Population propagation time\n";
    //cout<<transportstep<<"\n";

    cout<<"# End ########################################\n";
    excitonictransformed = 1;

}




void communicator_3rd::AddPattern(string istr)
{
//    this function setsup
//    coherentGSBK1  = icoherentGSBK1;
//    transportGSBK1 = itransportGSBK1;
//    coherentESEK1  = icoherentESEK1;
//    transportESEK1 = itransportESEK1;
//    coherentESAK1  = icoherentESAK1;
//    transportESAK1 = itransportESAK1;
//    coherentGSBK2  = icoherentGSBK2;
//    transportGSBK2 = itransportGSBK2;
//    coherentESEK2  = icoherentESEK2;
//    transportESEK2 = itransportESEK2;
//    coherentESAK2  = icoherentESAK2;
//    transportESAK2 = itransportESAK2;
//    coherentES1K3  = icoherentES1K3;
//    coherentES2K3  = icoherentES2K3;

    toolsIO tio;
    int numw;
    string*  singleitems = tio.toWords(numw, istr);

    if(numw == 1)
    {
        const char* istrc_str=istr.c_str();
        //parse letters
        coherentGSBK1  =  (istrc_str[0] == '1' ? 1 : 0);
        transportGSBK1 =  (istrc_str[1] == '1' ? 1 : 0);
        coherentESEK1  =  (istrc_str[2] == '1' ? 1 : 0);
        transportESEK1 =  (istrc_str[3] == '1' ? 1 : 0);
        coherentESAK1  =  (istrc_str[4] == '1' ? 1 : 0);
        transportESAK1 =  (istrc_str[5] == '1' ? 1 : 0);
        coherentGSBK2  =  (istrc_str[6] == '1' ? 1 : 0);
        transportGSBK2 =  (istrc_str[7] == '1' ? 1 : 0);
        coherentESEK2  =  (istrc_str[8] == '1' ? 1 : 0);
        transportESEK2 =  (istrc_str[9] == '1' ? 1 : 0);
        coherentESAK2  =  (istrc_str[10] == '1' ? 1 : 0);
        transportESAK2 =  (istrc_str[11] == '1' ? 1 : 0);
        coherentES1K3  =  (istrc_str[12] == '1' ? 1 : 0);
        coherentES2K3  =  (istrc_str[13] == '1' ? 1 : 0);
    }
    else
    {
        // parse words
            coherentGSBK1  =  (singleitems[0] == "1" ? 1 : 0);
            transportGSBK1 =  (singleitems[1] == "1" ? 1 : 0);
            coherentESEK1  =  (singleitems[2] == "1" ? 1 : 0);
            transportESEK1 =  (singleitems[3] == "1" ? 1 : 0);
            coherentESAK1  =  (singleitems[4] == "1" ? 1 : 0);
            transportESAK1 =  (singleitems[5] == "1" ? 1 : 0);
            coherentGSBK2  =  (singleitems[6] == "1" ? 1 : 0);
            transportGSBK2 =  (singleitems[7] == "1" ? 1 : 0);
            coherentESEK2  =  (singleitems[8] == "1" ? 1 : 0);
            transportESEK2 =  (singleitems[9] == "1" ? 1 : 0);
            coherentESAK2  =  (singleitems[10] == "1" ? 1 : 0);
            transportESAK2 =  (singleitems[11] == "1" ? 1 : 0);
            coherentES1K3  =  (singleitems[12] == "1" ? 1 : 0);
            coherentES2K3  =  (singleitems[13] == "1" ? 1 : 0);

    }


}



void communicator_3rd::SwitchConfigurations(
                                          int icoherentGSBK1,
                                          int itransportGSBK1,
                                          int icoherentESEK1,
                                          int itransportESEK1,
                                          int icoherentESAK1,
                                          int itransportESAK1,
                                          int icoherentGSBK2,
                                          int itransportGSBK2,
                                          int icoherentESEK2,
                                          int itransportESEK2,
                                          int icoherentESAK2,
                                          int itransportESAK2,
                                          int icoherentES1K3,
                                          int icoherentES2K3  )
{
    coherentGSBK1  = icoherentGSBK1;
    transportGSBK1 = itransportGSBK1;
    coherentESEK1  = icoherentESEK1;
    transportESEK1 = itransportESEK1;
    coherentESAK1  = icoherentESAK1;
    transportESAK1 = itransportESAK1;
    coherentGSBK2  = icoherentGSBK2;
    transportGSBK2 = itransportGSBK2;
    coherentESEK2  = icoherentESEK2;
    transportESEK2 = itransportESEK2;
    coherentESAK2  = icoherentESAK2;
    transportESAK2 = itransportESAK2;
    coherentES1K3  = icoherentES1K3;
    coherentES2K3  = icoherentES2K3;

}
void communicator_3rd::SwitchKI()
{
    coherentGSBK1  = 1;
    transportGSBK1 = 1;
    coherentESEK1  = 1;
    transportESEK1 = 1;
    coherentESAK1  = 1;
    transportESAK1 = 1;
}
void communicator_3rd::SwitchKII()
{
    coherentGSBK2  = 1;
    transportGSBK2 = 1;
    coherentESEK2  = 1;
    transportESEK2 = 1;
    coherentESAK2  = 1;
    transportESAK2 = 1;
}
void communicator_3rd::SwitchKIII()
{
    coherentES1K3  = 1;
    coherentES2K3  = 1;
}



void communicator_3rd::AddFluctuations(storage<asymptoticLF_complexv> imfun,storage<double> ikij,storage<double> igij)
{
	mfun = imfun;
	kij = ikij;
	gij = igij;
}
void communicator_3rd::AddFluctuations(storage<asymptoticLF_complexv> igfun,storage<asymptoticLF_complexv> imfun,storage<double> ikij,storage<double> igij)
{
	gfun = igfun;
	mfun = imfun;
	kij = ikij;
	gij = igij;
}



void communicator_3rd::SetEFields(dvector3d& ivE1, dvector3d& ivE2)
{
    vecE1 = ivE1;
    vecE2 = ivE2;
}
void communicator_3rd::SetEFields(dvector3d& ivE1, dvector3d& ivE2,dvector3d& ivE3,dvector3d& ivE4)
{
    vecE1 = ivE1;
    vecE2 = ivE2;
    vecE3 = ivE3;
    vecE4 = ivE4;
}
void communicator_3rd::SetAveraging(int& avera)
{
    averagingtype = avera;
}
void communicator_3rd::SetTemperature(double& tmpr)
{
    tempr = tmpr;
}


void communicator_3rd::SetupManifolds(int inumG, int inumE)
{
    numG = inumG;
    numE = inumE;
    numF = 0;
    numT = numG+numE+numF;
    ready = 1;
    band = 1;
    
}
void communicator_3rd::SetupManifolds(int inumG, int inumE, int inumF)
{
    numG = inumG;
    numE = inumE;
    numF = inumF;
    numT = numG+numE+numF;
    ready = 1;
    band = 1;
    
}



//------------------------------------

void communicator_3rd::Setup()
{
    // general flag
    ready = 0;
    excitonictransformed = 0;
    complexity =  "";
    band =1; 
    bath_complexity = "";
    experiment_complexity = "";


    // system parameters
    eigensys = true;
    evals.SetDimension(1);
    edips.SetDimension(2);
    dephasings.SetDimension(2);
    grPops.SetDimension(1); // ground manifold populations
    reorganizations.SetDimension(2);
    evec0.SetDimension(2);
    evec1.SetDimension(2);
    evec2.SetDimension(3);
	eten.SetDimension(2);
	ed_m.SetDimension(2);


    bosonic = 0;
    ham.SetDimension(2);
    anh.SetDimension(2);
    dip.SetDimension(1);
    dip2Corrections.SetDimension(1);
    coupl.SetDimension(2); // amplitudes for all chromophores [osc][chromophore]
	pos.SetDimension(1);
	d_m.SetDimension(1);


    // bath oscillators
    tempr = 0;
    spdf.SetDimension(1);
    mfun.SetDimension(1);
    gfun.SetDimension(1);
    cfun.SetDimension(1);

    gij.SetDimension(3); // amplitudes for all pairs
    kij.SetDimension(3); // amplitudes for all pairs
    zijkl.SetDimension(5);// correlation amplitudes for all hamiltonian elements
    numG = 0;
    numE = 0;
    numF = 0;
    numT = 0;
    nstates = 0;

    // transport rates:
    ratesg.SetDimension(2);
    ratese.SetDimension(2);
    nonsecular = 0;

    transport = 0; // turned off so far
    transportstep = 1.0;

    outputstring =  "";

    coherentGSBK1 = 1;
    transportGSBK1 = 1;
    coherentESEK1 = 1;
    transportESEK1 = 1;
    coherentESAK1 = 1;
    transportESAK1 = 1;
    coherentGSBK2 = 0;
    transportGSBK2 = 0;
    coherentESEK2 = 0;
    transportESEK2 = 0;
    coherentESAK2 = 0;
    transportESAK2 = 0;
    coherentES1K3 = 0;
    coherentES2K3 = 0;

    simType = 0;
    interactionpattern = "";

    vecE1 = dvector3d(0.,0.,1.);
    vecE2 = dvector3d(0.,0.,1.);
    vecE3 = dvector3d(0.,0.,1.);
    vecE4 = dvector3d(0.,0.,1.);
    averagingtype = 0;

	veck1=0.;
	veck2=0.;
	veck3=0.;
	veck4=0.;

	// more advanced parameters of the excitations (excitation central frequencies - in wavenumbers cm-1)
	xomega1=0.;
	xomega2=0.;
	xomega3=0.;
	xomega4=0.;

	// more advanced parameters of the excitations (excitation bandwidths - sigmas of Gaussians in wavenumbers cm-1)
	xsigma1=1.;
	xsigma2=1.;
	xsigma3=1.;
	xsigma4=1.;
    
    naturallinewidth = 0.;

	ham2sorter_l.SetDimension(1);
	ham2sorter_r.SetDimension(1);
}

void communicator_3rd::ReadFile(ifstream& ifs)
{
    // reading the  input file
    if(!ifs.is_open())
        cout<<"Error: communicator_3rd file has not been openned\n";
    else
    {
        ReadSystem(ifs);

        ReadBath(ifs);

        ReadExperiment(ifs);

        ReadSpecific(ifs);
        
        cout<<"all is read\n";

    }


}

void communicator_3rd::ReadSpecific(ifstream& ifs)
{
    cout<<"specific from communicator\n";
    return;
}

//void communicator_3rd::o_ReadSpecific(Open& input)
//{
//    return;
//}


void communicator_3rd::ReadFile(string& ifs)
{
    ifstream istr(ifs.c_str());
    if(istr.is_open())
    {
        ReadFile(istr);
    }
    else
    {
        cout<<"Error: input file cannot be read\n";
    }
}


void communicator_3rd::AddDephasings()
{
    int numt = numG+numE+numF;
    calculator_redfield calcR(evals);
    calcR.AddFluctuations(mfun,kij,gij);
    // next I make dephasings
    if(!dephasings.IsAlloc())
    {
        dephasings.Allocate(numt,numt);
        //cout<<"Dephasing rate matrix:\n";
        dephasings = calcR.AddLifetimeDephasings();
    }
//    if(!reorganizations.IsAlloc())
//    {
//        reorganizations.Allocate(numt,numt);
//        reorganizations = calcR.GetReorganizations();
//    }
    
}
void communicator_3rd::AddDephasingsAll()
{
    int numt = numG+numE+numF;
    calculator_redfield calcR(evals);
    calcR.AddFluctuations(mfun,kij,gij);
    // next I make dephasings
    if(!dephasings.IsAlloc())
    {
        dephasings.Allocate(numt,numt);
        //cout<<"Dephasing rate matrix:\n";
        calcR.AddLifetimeDephasings();
        dephasings = calcR.AddPureDephasings();
    }
//    if(!reorganizations.IsAlloc())
//    {
//        reorganizations.Allocate(numt,numt);
//        reorganizations = calcR.GetReorganizations();
//    }
    
}

void communicator_3rd::AddTransport(storage<double>& iratesG,storage<double>& iratesE,double tstep)
{
    ratesg = iratesG;
    ratese = iratesE;
    transport = 1;
    transportstep = tstep;
}



void communicator_3rd::AddTransport()
{
    
    cout<<"Error: do not call function \"communicator_3rd::AddTransport()\"!!!\n";
    
    // here I create redfield rates and
    // dephasing constants from mfun

    transport = 1;
    
    storage<asymptoticLF_complexv> g1un; // gfunctions first derivatives
    
    std::size_t found;
    found = experiment_complexity.find("ModRed");
    if(found!=std::string::npos)
	{
        int numBosc = spdf.GetSize();
        if(!g1un.IsAlloc())
            g1un.Allocate(numBosc);
        for(int ind=0; ind<numBosc; ind++)
        {
            numericalSD spd(tempr,spdf.data1D[ind]);
            g1un.data1D[ind] = spd.GetGfD1f();
        }
        
    }

    // GG block
    if(true)
    {
        ratesg.Allocate(numG,numG);

        // block-specific data transform
        cout<<"Redfield rate matrix: GG:\n";
        int numo = mfun.GetSize();
        storage<double> ens(1);
        storage<double> ekij(3);
        storage<double> egij(3);
        ens.Allocate(numG);
        ekij.Allocate(numG,numG,numo);
        egij.Allocate(numG,numG,numo);
        for(int is=0; is<numG; is++)
            ens.data1D[is]=evals.data1D[is];
        for(int io=0; io<numo; io++)
        for(int is=0; is<numG; is++)
        for(int iz=0; iz<numG; iz++)
        {
            ekij.data3D[is][iz][io]=kij.data3D[is][iz][io];
            egij.data3D[is][iz][io]=gij.data3D[is][iz][io];
        }
        calculator_redfield calcR(ens);
        calcR.AddFluctuations(mfun,ekij,egij);
        calcR.AddMijkl(zijkl);

        
        found = experiment_complexity.find("ModRed");
        if(found!=std::string::npos)
        {
            calcR.gdun = g1un;
            ratesg = calcR.GetTransportRatesModRed();
        }
        else
            ratesg = calcR.AddTransportRates();
    }
    // EE block
    if(true)
    {
        ratese.Allocate(numE,numE);

        // block-specific data transform
        cout<<"Redfield rate matrix: EE:\n";
        int numo = mfun.GetSize();
        storage<double> ens(1);
        storage<double> ekij(3);
        storage<double> egij(3);
        ens.Allocate(numE);
        ekij.Allocate(numE,numE,numo);
        egij.Allocate(numE,numE,numo);
        for(int is=0; is<numE; is++)
            ens.data1D[is]=evals.data1D[is+numG];
        for(int io=0; io<numo; io++)
            for(int is=0; is<numE; is++)
                for(int iz=0; iz<numE; iz++)
                {
                    ekij.data3D[is][iz][io]=kij.data3D[is+numG][iz+numG][io];
                    egij.data3D[is][iz][io]=gij.data3D[is+numG][iz+numG][io];
                }
        calculator_redfield calcR(ens);
        calcR.AddFluctuations(mfun,ekij,egij);
        calcR.AddMijkl(zijkl);
        
        found = experiment_complexity.find("ModRed");
        if(found!=std::string::npos)
        {
            calcR.gdun = g1un;
            ratese = calcR.GetTransportRatesModRed();
        }
        else
            ratese = calcR.AddTransportRates();
    }
    // these are rates from right to left

}

void communicator_3rd::ReadExperiment(ifstream& ifs)
{
    if(band == 1)
        ReadExperimentLin(ifs);
    else
        ReadExperiment3rd(ifs);
}

void communicator_3rd::ReadExperimentLin(ifstream& ifs)
{
    storage<double> arr(2); // two dimensional
    storage<int> iri(2); // two dimensional

    toolsIO tio;
    string str = "";
//                 cout<<"Reading experiment parameters\n";
   
    cout<<"Reading experiment setup line:\n";
    str = "";
    tio.StreamSkipTrailers(&ifs);
    getline(ifs,str);
    tio.StreamSkipTrailers(str);
    experiment_complexity = str;
    cout<<"Got: \""<<str<<"\"\n";

    averagingtype = 1;
    std::size_t found;
    found = experiment_complexity.find("Isotropic");
    if(found!=std::string::npos)
       averagingtype = 0;

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentPolarizations:",
		"Reading experiment laser polarizations\n", 
		2,3, arr,inputfile.str()))
	{
        vecE1 = dvector3d(arr.data2D[0][0],arr.data2D[0][1],arr.data2D[0][2]);
        vecE2 = dvector3d(arr.data2D[1][0],arr.data2D[1][1],arr.data2D[1][2]);
	}
    
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentWavevectors:",
		"Reading experiment laser propagation directions\n", 
		2,3, arr,inputfile.str()))
	{
        veck1 = dvector3d(arr.data2D[0][0],arr.data2D[0][1],arr.data2D[0][2]);
        veck2 = dvector3d(arr.data2D[1][0],arr.data2D[1][1],arr.data2D[1][2]);
	}
    
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInteractionFrequencies:",
		"Reading experiment laser interaction frequencies\n", 
		2,1, arr,inputfile.str()))
	{
        xomega1 = arr.data2D[0][0];
        xomega2 = arr.data2D[1][0];
	}
    //cout<<"reading experiment parameters\n";
    
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInitialFrequency:",
		"Reading initial FFT frequency\n", 
		1,1, arr,inputfile.str()))
	{
        ifre1 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentFinalFrequency:",
		"Reading final FFT frequency\n", 
		1,1, arr,inputfile.str()))
	{
        ffre1 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<int>(
		"ExperimentNumberOfPoints:",
		"Reading the number of points in the interval\n", 
		1,1, iri,inputfile.str()))
	{
        nump = iri.data2D[0][0];
	}

        
        
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentNaturalLinewidth:",
		"Reading the natural linewidth parameter\n", 
		1,1, arr,inputfile.str()))
	{
        naturallinewidth = arr.data2D[0][0];
	}

        
}


void communicator_3rd::ReadExperiment3rd(ifstream& ifs)
{
                toolsIO tio;
                string str = "";

	storage<double> arr(2);
	storage<int> iri(2);
            std::size_t found;


    cout<<"Reading experiment setup line:\n";
    str = "";
    tio.StreamSkipTrailers(&ifs);
    getline(ifs,str);
    tio.StreamSkipTrailers(str);
    experiment_complexity = str;
    cout<<"Got: \""<<str<<"\"\n";


        found = experiment_complexity.find("2DES");
        if(found!=std::string::npos)
	{
		// very nice; correct 
		;
        }
	else 
	{
		cout<<"Error: wrong experiment - 2DES?\n";
	}


	simType=0;
        found = experiment_complexity.find("WTW");
        if(found!=std::string::npos)
	{
		simType= 2;
        }
        found = experiment_complexity.find("TWW");
        if(found!=std::string::npos)
	{
		simType= 1;
        }
        found = experiment_complexity.find("WWT");
        if(found!=std::string::npos)
	{
		simType= 3;
        }
        found = experiment_complexity.find("TTT");
        if(found!=std::string::npos)
	{
		simType= 4;
        }
	if(simType == 0)
		std::cout<<"Error in communicator_3rd: simulation type is missing. Further behavior is unpredictable\n";



    // swithing on and off various diagrams
    // the code is as follows:
    // CGSBK1 TGSBK1 CESEK1 TESEK1 CESAK1 TESAK1 CGSBK2 TGSBK2 CESEK2 TESEK2 CESAK2 TESAK2 CES1K3 CES2K3
	if(tio.LookUpAndReadSquareMatrix<int>(
		"ExperimentDiagramTypes:",
		"Reading experiment diagram patterm\n", 
		1,14, iri,inputfile.str()))
	{
		SwitchConfigurations(iri.data2D[0][0],iri.data2D[0][1],
		iri.data2D[0][2],iri.data2D[0][3],
		iri.data2D[0][4],iri.data2D[0][5],
		iri.data2D[0][6],iri.data2D[0][7],
		iri.data2D[0][8],iri.data2D[0][9],
		iri.data2D[0][10],iri.data2D[0][11],
		iri.data2D[0][12],iri.data2D[0][13]);
	}

    averagingtype = 1;
    found = experiment_complexity.find("Isotropic");
    if(found!=std::string::npos)
       averagingtype = 0;


	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentPolarizations:",
		"Reading experiment laser polarizations\n", 
		4,3, arr,inputfile.str()))
	{
        vecE1 = dvector3d(arr.data2D[0][0],arr.data2D[0][1],arr.data2D[0][2]);
        vecE2 = dvector3d(arr.data2D[1][0],arr.data2D[1][1],arr.data2D[1][2]);
        vecE3 = dvector3d(arr.data2D[2][0],arr.data2D[2][1],arr.data2D[2][2]);
        vecE4 = dvector3d(arr.data2D[3][0],arr.data2D[3][1],arr.data2D[3][2]);
	}
    
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentWavevectors:",
		"Reading experiment laser propagation directions\n", 
		4,3, arr,inputfile.str()))
	{
        veck1 = dvector3d(arr.data2D[0][0],arr.data2D[0][1],arr.data2D[0][2]);
        veck2 = dvector3d(arr.data2D[1][0],arr.data2D[1][1],arr.data2D[1][2]);
        veck3 = dvector3d(arr.data2D[2][0],arr.data2D[2][1],arr.data2D[2][2]);
        veck4 = dvector3d(arr.data2D[3][0],arr.data2D[3][1],arr.data2D[3][2]);
	}
    
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInteractionFrequencies:",
		"Reading experiment laser interaction frequencies\n", 
		4,1, arr,inputfile.str()))
	{
        xomega1 = arr.data2D[0][0];
        xomega2 = arr.data2D[1][0];
        xomega3 = arr.data2D[2][0];
        xomega4 = arr.data2D[3][0];
	}


	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInitialFrequency1:",
		"Reading initial FFT frequency 1\n", 
		1,1, arr,inputfile.str()))
	{
        ifre1 = arr.data2D[0][0];
	}

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentFinalFrequency1:",
		"Reading final FFT frequency 1\n", 
		1,1, arr,inputfile.str()))
	{
        ffre1 = arr.data2D[0][0];
	}

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInitialFrequency2:",
		"Reading initial FFT frequency 2\n", 
		1,1, arr,inputfile.str()))
	{
        ifre2 = arr.data2D[0][0];
	}

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentFinalFrequency2:",
		"Reading final FFT frequency 2\n", 
		1,1, arr,inputfile.str()))
	{
        ffre2 = arr.data2D[0][0];
	}

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInitialFrequency3:",
		"Reading initial FFT frequency 3\n", 
		1,1, arr,inputfile.str()))
	{
        ifre3 = arr.data2D[0][0];
	}

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentFinalFrequency3:",
		"Reading final FFT frequency 3\n", 
		1,1, arr,inputfile.str()))
	{
        ffre3 = arr.data2D[0][0];
	}

	/////////////////////////
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentDelay3:",
		"Reading time delay 3\n", 
		1,1, arr,inputfile.str()))
	{
        tf3 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentFinalDelay3:",
		"Reading final time delay 3\n", 
		1,1, arr,inputfile.str()))
	{
        tf3 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInitialDelay3:",
		"Reading initial time delay 3\n", 
		1,1, arr,inputfile.str()))
	{
        ti3 = arr.data2D[0][0];
	}

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentDelay2:",
		"Reading time delay 2\n", 
		1,1, arr,inputfile.str()))
	{
        tf2 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentFinalDelay2:",
		"Reading final time delay 2\n", 
		1,1, arr,inputfile.str()))
	{
        tf2 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInitialDelay2:",
		"Reading initial time delay 2\n", 
		1,1, arr,inputfile.str()))
	{
        ti2 = arr.data2D[0][0];
	}

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentDelay1:",
		"Reading time delay 1\n", 
		1,1, arr,inputfile.str()))
	{
        tf1 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentFinalDelay1:",
		"Reading final time delay 1\n", 
		1,1, arr,inputfile.str()))
	{
        tf1 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentInitialDelay1:",
		"Reading initial time delay 1\n", 
		1,1, arr,inputfile.str()))
	{
        ti1 = arr.data2D[0][0];
	}

	////////////////////////////
	if(tio.LookUpAndReadSquareMatrix<int>(
		"ExperimentNumberOfPoints:",
		"Reading the number of points in the interval\n", 
		1,1, iri,inputfile.str()))
	{
        nump = iri.data2D[0][0];
	}

	if(simType == 4)
	{

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentResonantFrequency1:",
		"Reading experiment resonant frequency 1\n", 
		1,1, arr,inputfile.str()))
	{
        ifre1 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentResonantFrequency2:",
		"Reading experiment resonant frequency 2\n", 
		1,1, arr,inputfile.str()))
	{
        ifre2 = arr.data2D[0][0];
	}
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentResonantFrequency3:",
		"Reading experiment resonant frequency 3\n", 
		1,1, arr,inputfile.str()))
	{
        ifre3 = arr.data2D[0][0];
	}
	}
	
	
	
	
	
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExperimentNaturalLinewidth:",
		"Reading the natural linewidth parameter\n", 
		1,1, arr,inputfile.str()))
	{
        naturallinewidth = arr.data2D[0][0];
	}


    
}




void communicator_3rd::ReadSystem(ifstream& ifs)
{
    toolsIO tio;
    string str = "";


    cout<<"Reading setup line:\n";
    str = "";
    tio.StreamSkipTrailers(&ifs);
    getline(ifs,str);
    tio.StreamSkipTrailers(str);
    complexity = str;
    cout<<"Got: \""<<complexity<<"\"\n";

    std::size_t found;
    found = complexity.find("Exciton");
    if(found!=std::string::npos)
	{
        eigensys = false;

        bosonic = false;
        std::size_t found = complexity.find("Bosonic");
        if(found!=std::string::npos)
            bosonic = true;
	}
    else
        eigensys = true;

    found = complexity.find("2band");
    if(found!=std::string::npos)
        band = 2;
    else
        band = 1;

    //found = complexity.find("nonsecular");
    //if(found!=std::string::npos)
    //    nonsecular = 1;


	// reading new form
	storage<double> arr(2);
	storage<int> iri(2);



	// eigensys part
	if(band == 2)
	if(tio.LookUpAndReadSquareMatrix<int>(
		"EigensysNumberOf012Levels:",
		"Reading the number of zero, 1 and 2-band levels\n", 
		3,1, iri,inputfile.str()))
	{
                numG = iri.data2D[0][0];
                numE = iri.data2D[1][0];
                numF = iri.data2D[2][0];
		numT = numG+numE+numF;
	}

	if(band == 1)
	if(tio.LookUpAndReadSquareMatrix<int>(
		"EigensysNumberOf01Levels:",
		"Reading the number of zero and 1-band levels\n", 
		2,1, iri,inputfile.str()))
	{
                numG = iri.data2D[0][0];
                numE = iri.data2D[1][0];
                numF = 0;
		numT = numG+numE+numF;
	}


	if(tio.LookUpAndReadSquareMatrix<double>(
		"EigensysLevelEnergies:",
		"Reading the energies of all energy levels\n", 
		numT,1, arr,inputfile.str()))
	{
        	evals.Allocate(numT);
		for(int ind=0; ind<numT; ind ++)
			evals.data1D[ind]=arr.data2D[ind][0];
	}


	if(tio.LookUpAndReadSquareMatrix<int>(
		"EigensysTransitionDipoles:",
		"Reading the transition dipoles between all energy levels\n", 
		1,1, iri,inputfile.str()))
	{	// this one is for the number of entries
		int numX = iri.data2D[0][0];

		// next the all entries

		string sstr = inputfile.str();
		string str = "";
		// looking for entry
		std::size_t found = sstr.find("EigensysTransitionDipoles:");
		sstr = sstr.substr(found);
		std::istringstream ifs(sstr);
		getline(ifs,str);
		getline(ifs,str);
		if(arr.IsAlloc()) arr.Delete();
        	arr.Allocate(numX,5);
		if(tio.ReadRectangular(arr, ifs))
                {
                      cout<<"Error: Some problem with file reading\n";
                }		// assignment
        	edips.Allocate(numT,numT);
		int i2,i1;
		for(int ind=0; ind<numX; ind ++){
			i1 = (int)arr.data2D[ind][0];
			i2 = (int)arr.data2D[ind][1];
			edips.data2D[i1][i2]=dvector3d(arr.data2D[ind][2],arr.data2D[ind][3],arr.data2D[ind][4]);
			edips.data2D[i2][i1]=edips.data2D[i1][i2];
		}
	}

	if(tio.LookUpAndReadSquareMatrix<int>(
		"EigensysMagneticDipoles:",
		"Reading the magnetic transition dipoles between all energy levels\n", 
		1,1, iri,inputfile.str()))
	{	// this one is for the number of entries
		int numX = iri.data2D[0][0];

		// next the all entries

		string sstr = inputfile.str();
		string str = "";
		// looking for entry
		std::size_t found = sstr.find("EigensysMagneticDipoles:");
		sstr = sstr.substr(found);
		std::istringstream ifs(sstr);
		getline(ifs,str);
		getline(ifs,str);
		if(arr.IsAlloc()) arr.Delete();
        	arr.Allocate(numX,5);
		if(tio.ReadRectangular(arr, ifs))
                {
                      cout<<"Error: Some problem with file reading\n";
                }		// assignment
        	ed_m.Allocate(numT,numT);
		int i2,i1;
		for(int ind=0; ind<numX; ind ++){
			i1 = (int)arr.data2D[ind][0];
			i2 = (int)arr.data2D[ind][1];
			ed_m.data2D[i1][i2]=dvector3d(arr.data2D[ind][2],arr.data2D[ind][3],arr.data2D[ind][4]);
			ed_m.data2D[i2][i1]=ed_m.data2D[i1][i2];
		}
	}


	if(tio.LookUpAndReadSquareMatrix<int>(
		"EigensysQuadrupoles:",
		"Reading the magnetic transition dipoles between all energy levels\n", 
		1,1, iri,inputfile.str()))
	{	// this one is for the number of entries
		int numX = iri.data2D[0][0];

		// next the all entries

		string sstr = inputfile.str();
		string str = "";
		// looking for entry
		std::size_t found = sstr.find("EigensysQuadrupoles:");
		sstr = sstr.substr(found);
		std::istringstream ifs(sstr);
		getline(ifs,str);
		getline(ifs,str);
		if(arr.IsAlloc()) arr.Delete();
        	arr.Allocate(numX,8);
		if(tio.ReadRectangular(arr, ifs))
                {
                      cout<<"Error: Some problem with file reading\n";
                }		// assignment
        	eten.Allocate(numT,numT);
		int i2,i1;
		for(int ind=0; ind<numX; ind ++){
			i1 = (int)arr.data2D[ind][0];
			i2 = (int)arr.data2D[ind][1];
			eten.data2D[i1][i2]=dtensor3x3(
                arr.data2D[ind][2],//xx
                arr.data2D[ind][3],//xy
                arr.data2D[ind][4],//xz
                arr.data2D[ind][3],//yx
                arr.data2D[ind][5],//yy
                arr.data2D[ind][6],//yz
                arr.data2D[ind][4],//zx
                arr.data2D[ind][6],//zy
                arr.data2D[ind][7] //zz                
                                          );
			eten.data2D[i2][i1]=eten.data2D[i1][i2];
		}
	}



	// excitonic part

	// excitonic number of sites 
	if(tio.LookUpAndReadSquareMatrix<int>(
		"ExcNumberOfSites:",
		"Reading excitonic number of sites\n", 
		1, 1, iri,inputfile.str()))
	{
		numE = iri.data2D[0][0];
		numG = 1;
	}

	// excitonic main hamiltonian
	if(tio.LookUpAndReadTriangularMatrix<double>(
		"ExcMainHamiltonian:",
		"Reading main excitonic hamiltonian\n", 
		numE, arr,inputfile.str()))
	{
		ham.Allocate(numE,numE);
                for(int inl=0; inl<numE; inl ++)
                for(int ind=0; ind<=inl; ind ++){
			ham.data2D[inl][ind] = arr.data2D[inl][ind];
			ham.data2D[ind][inl] = ham.data2D[inl][ind] ;
		}
	}


	// excitonic biparticle anharmonicities
    if(band == 2)
	if(tio.LookUpAndReadTriangularMatrix<double>(
		"ExcBiparticleAnharmonicites:",
		"Reading excitonic anharmonicities\n", 
		numE, arr,inputfile.str()))
	{
		anh.Allocate(numE,numE);
                for(int inl=0; inl<numE; inl ++)
                for(int ind=0; ind<=inl; ind ++){
			anh.data2D[inl][ind] = arr.data2D[inl][ind];
			anh.data2D[ind][inl] = anh.data2D[inl][ind] ;
		}
	}


	// excitonic electric vectors
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExcElectricDipoles:",
		"Reading excitonic electric dipoles\n", 
		numE, 3, arr,inputfile.str()))
	{
		dip.Allocate(numE);
                // making dipoles
                for(int ind=0; ind<numE; ind ++)
                      dip.data1D[ind]=dvector3d(arr.data2D[ind][0],arr.data2D[ind][1],arr.data2D[ind][2]);
	}

	// excitonic position vectors
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExcPositionVectors:",
		"Reading excitonic position vectors\n", 
		numE, 3, arr,inputfile.str()))
	{
		pos.Allocate(numE);
                // making dipoles
                for(int ind=0; ind<numE; ind ++)
                      pos.data1D[ind]=dvector3d(arr.data2D[ind][0],arr.data2D[ind][1],arr.data2D[ind][2]);
	}

	// reading excitonic magnetic dipoles
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExcMagneticDipoles:",
		"Reading excitonic magnetic dipoles\n", 
		numE, 3, arr,inputfile.str()))
	{
		d_m.Allocate(numE);
                // making dipoles
                for(int ind=0; ind<numE; ind ++)
                      d_m.data1D[ind]=dvector3d(arr.data2D[ind][0],arr.data2D[ind][1],arr.data2D[ind][2]);
	}


	// excitonic anharmonic vectors
    if(band == 2)
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExcElectricDipoleAnharmonicities:",
		"Reading excitonic electric dipole anharmonicities\n", 
		numE, 3, arr,inputfile.str()))
	{
		dip2Corrections.Allocate(numE);
                // making dipoles
                for(int ind=0; ind<numE; ind ++)
                      dip2Corrections.data1D[ind]=dvector3d(arr.data2D[ind][0],arr.data2D[ind][1],arr.data2D[ind][2]);
	}



	// finalizing setup
	if(!eigensys){
		numG = 1;
		numF = 0; 
		if(band == 2){
			if(bosonic)
				numF = numE*(numE+1)/2;
			else
				numF = numE*(numE-1)/2;
		}
		numT = numG+numE+numF;
        }

}

void communicator_3rd::ReadBath(ifstream& ifs)
{
    string str = "";

    toolsIO tio;

    int numBosc;
	int numl,numr;

	storage<double> arr(2);
    storage<int> iri(2);


    cout<<"Reading bath setup line:\n";
    str = "";
    tio.StreamSkipTrailers(&ifs);
    getline(ifs,str);
    tio.StreamSkipTrailers(str);
    bath_complexity = str;
    cout<<"Got: \""<<bath_complexity<<"\"\n";




	// bath temperature
	if(tio.LookUpAndReadSquareMatrix<double>(
		"BathTemperature:",
		"Reading bath temperature\n", 
		1, 1, arr,inputfile.str()))
	{
		tempr = arr.data2D[0][0];
	}

	// bath oscillator number
	if(tio.LookUpAndReadSquareMatrix<int>(
		"BathNumberOfOscillators:",
		"Reading number of bath oscillators\n", 
		1, 1, iri,inputfile.str()))
	{
		numBosc = iri.data2D[0][0];
	}

	if(!LookUpAndReadBathSpectraldensities(
		"BathSpectralDensities:",
		"Reading all bath spectral densities\n", numBosc))
	{
		std::cout<<"Warning: bath spectral densities are missing: the code may not finish correctly\n";
	}


    if(eigensys)
    {
	string sstr = inputfile.str();
	string str;


	// looking for entry
	std::size_t found = sstr.find("BathSystemCorrelatedCouplingMagnitudes:");
	if(found != std::string::npos)
	{
		kij.Allocate(numT,numT,numBosc);
		gij.Allocate(numT,numT,numBosc);

       		sstr = sstr.substr(found);
		std::istringstream ifs(sstr);
		getline(ifs,str);
		cout<<"Reading bath-system coupling matrices\n";
		for(int inb = 0; inb<numBosc; inb++)
		{
			if(arr.IsAlloc()) arr.Delete();
	                arr.Allocate(numT,numT);
	                if(tio.ReadTriangular(arr, ifs))
			{
				cout<<"Error: Some problem with file reading\n";
				return;
			}
			for(int indl=0; indl<numT; indl++)
			for(int indr=0; indr<=indl; indr++)
			{
				kij.data3D[indr][indl][inb] = arr.data2D[indl][indr] ;
				kij.data3D[indl][indr][inb] = arr.data2D[indl][indr] ;
			}

			if(arr.IsAlloc()) arr.Delete();
	                arr.Allocate(numT,numT);
	                if(tio.ReadTriangular(arr, ifs))
			{
				cout<<"Error: Some problem with file reading\n";
				return;
			}
			for(int indl=0; indl<numT; indl++)
			for(int indr=0; indr<=indl; indr++)
			{
				gij.data3D[indr][indl][inb] = arr.data2D[indl][indr] ;
				gij.data3D[indl][indr][inb] = arr.data2D[indl][indr] ;
			}
		}
	}
    }
    else
    {
	// these are purely excitonic properties        
        //cout<<"reading couplings to the bath\n";
        
        int ctype=0;
        // reading the flag if fluctuations are correlated or not
        //if(complexity.compare("key-setup-complete")==0)
        {
            std::size_t found;
            ctype = 0;
            found = bath_complexity.find("Correlated");
            if(found!=std::string::npos)
                ctype = 1;
            found = bath_complexity.find("UnCorrelated");
            if(found!=std::string::npos)
                ctype = 0;

        }
        
        

/*
	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExcBathSystemCorrelatedCouplingConstants:",
		"Reading bath-system coupling constants. Coupling term is |n><n|*w*c*q. The reorganization energy of that mode is (wc^2)/2. Now the spectral density file contains these c numbers. The whole \"spectral density\" is scaled by this coupling parameter.\n", 
		numBosc,numE, arr,inputfile.str()))
	{
		coupl.Allocate(numBosc,numE);
                for(int inl=0; inl<numBosc; inl ++)
                for(int ind=0; ind<numE; ind ++)
			coupl.data2D[inl][ind] = arr.data2D[inl][ind];
	}
*/

        
        if(ctype == 0)
        {
           // uncorrelated 
		numl=numBosc;
		numr=numE;

	if(tio.LookUpAndReadSquareMatrix<double>(
		"ExcBathSystemCouplingMagnitudes:",
		"Reading bath-system coupling amplitudes\n", 
		numl,numr, arr,inputfile.str()))
	{
		coupl.Allocate(numl,numr);
                for(int inl=0; inl<numl; inl ++)
                for(int ind=0; ind<numr; ind ++)
			coupl.data2D[inl][ind] = arr.data2D[inl][ind];
	}


	}
        else
        {
		//correlated
            	numl = numE;
		numr = numE;
		storage<double> cmatrix(3);
		cmatrix.Allocate(numBosc,numE,numE);
		for(int inb = 0; inb<numBosc; inb++)
		{

			stringstream ss;
			ss<<inb;
        
//            str = "";
	if(tio.LookUpAndReadTriangularMatrix<double>(
		"ExcBathSystemCorrelatedCouplingMagnitudes_"+ss.str()+":",
		"Reading bath-system coupling amplitudes "+ss.str()+"\n", 
		numl, arr,inputfile.str()))
	{
		for(int inl=0; inl<numl; inl ++)
                for(int ind=0; ind<=inl; ind ++){
			cmatrix.data3D[inb][inl][ind] = arr.data2D[inl][ind];
			cmatrix.data3D[inb][ind][inl] = arr.data2D[inl][ind];
		}
	}

		}
            coupl = cmatrix;
            // deallocation is automatic
           
        }
       
    }


}


int communicator_3rd::LookUpAndReadBathSpectraldensities(string keyword,string label, int numBosc){

	string sstr = inputfile.str();
	string str = "";

	// looking for entry
	std::size_t found = sstr.find(keyword);
	if(found != std::string::npos)
	{
        	sstr = sstr.substr(found);
		std::istringstream ifs(sstr);
		getline(ifs,keyword);
           	toolsIO tio;
		cout<<label;



    
    
    spdf.Allocate(numBosc);
    
    mfun.Allocate(numBosc);
    
    gfun.Allocate(numBosc);

    cfun.Allocate(numBosc);
    
    
    for(int ind = 0; ind<numBosc; ind++)
    {
        int reading = 1;
        
        // taking filenames from the ifile
        str = "";
        tio.StreamSkipTrailers(ifs);
        getline(ifs,str);
        tio.StreamSkipTrailers(str);
        cout<<"Reading and transforming file: "<<str<<"\n";
        
        numericalSD spd(tempr);
//        spd.SetZPLGamma(gammaZPL);
        
        
        if(true)
        {
            // checking the wrk files
            string fstr = str + ".wrk";
            ifstream ifsl(fstr.c_str());
            if (ifsl.is_open())
            {
                // reading backup .wrk file
                spdf.data1D[ind].ReadF(&ifsl);
            }
            else
            {
                // reading text and making backup
                ifstream ifsl(str.c_str());
                interpolationF<double> readf;
                readf.ReadF(&ifsl);
                spdf.data1D[ind] = readf;
                spdf.data1D[ind].causality = 3;
                spdf.data1D[ind].SetTitle(str);
                ifsl.close();

                // making backup
                spdf.data1D[ind].SaveF(fstr);
            }
            spd.SetSD(spdf.data1D[ind]);
        }
        if(true)
        {
            // checking the wrk files
            string fstr = str + ".mfun.wrk";
            ifstream ifsl(fstr.c_str());
            if (ifsl.is_open())
            {
                // reading backup
                mfun.data1D[ind].ReadF(&ifsl);
            }
            else
            {
                cout<<"# Making M function\n";
                // making backup
                mfun.data1D[ind] = spd.GetMf();// the whole function is returned

                // making backup
                mfun.data1D[ind].SaveF(fstr);
            }
        }
        if(true)
        {
            // checking the wrk files
            string fstr = str + ".gfun.wrk";
            ifstream ifsl(fstr.c_str());
            if (ifsl.is_open())
            {
                // reading backup
                gfun.data1D[ind].ReadF(&ifsl);
            }
            else
            {
                cout<<"# Making G function\n";
                // making backup
                gfun.data1D[ind] = spd.GetGf();// the whole function is returned
                // making backup
                gfun.data1D[ind].SaveF(fstr);
            }
        }
        if(true)
        {
            // checking the wrk files
            string fstr = str + ".cfun.wrk";
            ifstream ifsl(fstr.c_str());
            if (ifsl.is_open())
            {
                // reading backup
                cfun.data1D[ind].ReadF(&ifsl);
            }
            else
            {
                cout<<"# Making C function\n";
                // making backup
                cfun.data1D[ind] = spd.GetCf();// the whole function is returned
                // making backup
                cfun.data1D[ind].SaveF(fstr);
            }
        }
        cout<<"Done processing file: "<<str<<"\n";
  }
        
        
    
	return 1;
	}
	return 0;
}


///////// publishers


void communicator_3rd::Publish(string& ife)
{
    ofstream ifs(ife.c_str());
    Publish(ifs);
    ifs.close();
}



void communicator_3rd::Publish(ofstream& ofs)
{
    
     cout<<"Publishing result\n";
    if(simType == 0) // for linear spectra
    {
        ofs<<outputstring;
        ofs.precision(12);
        double fre;
        double dfre = (ffre1-ifre1)/nump;
        for(int ind = 0; ind<nump; ind++)
        {
            fre = ifre1+dfre*ind;
            ofs<<fre<<"\t"<<fre*(signal1d.Get(fre)).imag()<<"\t"<<fre*(signal1d.Get(fre)).real()<<"\n";
            //if(ind == nump-1)
            //{
            //    cout.precision(16);
            //    cout<<fre<<" "<<signal.Get(fre)<<"\n";
            //}
        }
    }   

    if(simType == 1)
    {
        // this is only for 2Q2D
        ofs<<outputstring;
        ofs.precision(12);

        double fre2;
        double dfre2 = (ffre2-ifre2)/nump;
        double fre3;
        double dfre3 = (ffre3-ifre3)/nump;
        for(int ind1 = 0; ind1<nump; ind1++)
            for(int ind2 = 0; ind2<nump; ind2++)
            {
            fre2 = ifre2+dfre2*ind1;
            fre3 = ifre3+dfre3*ind2;
            complexv result = signal.Get(fre2,fre3);
            ofs<<fre3<<"\t"<<fre2<<"\t"<<result.real()<<"\t";
            ofs<<result.imag()<<"\n";
        }
    }
    else if(simType == 2)
    {
        // this is only for rephasing and non-rephasing 2D
        ofs<<outputstring;
        ofs.precision(12);

        double fre1;
        double dfre1 = (ffre1-ifre1)/nump;
        double fre3;
        double dfre3 = (ffre3-ifre3)/nump;
        for(int ind1 = 0; ind1<nump; ind1++)
            for(int ind3 = 0; ind3<nump; ind3++)
            {
            fre1 = ifre1+dfre1*ind1;
            fre3 = ifre3+dfre3*ind3;
            complexv result = signal.Get(fre1,fre3);
            ofs<<fre3<<"\t"<<fre1<<"\t"<<result.real()<<"\t";
            ofs<<result.imag()<<"\n";
        }
    }
    else if(simType == 3)
    {
        // this is only for 2Q2D
        ofs<<outputstring;
        ofs.precision(12);

        double fre1;
        double dfre1 = (ffre1-ifre1)/nump;
        double fre2;
        double dfre2 = (ffre2-ifre2)/nump;
        for(int ind1 = 0; ind1<nump; ind1++)
            for(int ind2 = 0; ind2<nump; ind2++)
            {
            fre1 = ifre1+dfre1*ind1;
            fre2 = ifre2+dfre2*ind2;
            complexv result = signal.Get(fre1,fre2);
            ofs<<fre2<<"\t"<<fre1<<"\t"<<result.real()<<"\t";
            ofs<<result.imag()<<"\n";
        }
    }
    else if(simType == 4)
    {
        PublishSpecial4(ofs);
    }
}

void communicator_3rd::Publish2DRe(ofstream& ofs)
{

    //ofs<<"# QCFP communicator_3rd_excitons class signal Re\n";
    ofs.precision(12);


    storage<double> tstr(2);
    tstr.Allocate(nump,nump);


    storage<complexv> result(2);
    result.data2D = signal.DirectExchange(result.data2D);

    for(int ind1 = 0; ind1<nump; ind1++)
    for(int ind3 = 0; ind3<nump; ind3++)
    {
        tstr.data2D[ind1][ind3] = result.data2D[ind1][ind3].real();
    }
    result.data2D = signal.DirectExchange(result.data2D);



    toolsIO tio;
    tio.WriteRectangular(&tstr, ofs);
}

void communicator_3rd::Publish2DIm(ofstream& ofs)
{

    //ofs<<"# QCFP communicator_3rd_excitons class signal Im\n";
    ofs.precision(12);


    storage<double> tstr(2);
    tstr.Allocate(nump,nump);


    storage<complexv> result(2);
    result.data2D = signal.DirectExchange(result.data2D);

    for(int ind1 = 0; ind1<nump; ind1++)
    for(int ind3 = 0; ind3<nump; ind3++)
    {
        tstr.data2D[ind1][ind3] = result.data2D[ind1][ind3].imag();
    }
    result.data2D = signal.DirectExchange(result.data2D);



    toolsIO tio;
    tio.WriteRectangular(&tstr, ofs);
}


void communicator_3rd::FinalizeESSystem(int nonmarkovian)
{
    int numBosc = mfun.GetSize();
    
        // shift all energies because of the coupling to the bath
        for(int ind=0;ind<numT;ind++){

            double reor = 0.0;
            for(int ib=0;ib<numBosc;ib++)
                reor += -(gij.data3D[ind][ind][ib]*mfun.data1D[ib].Get(0.0)).imag();
                
            evals.data1D[ind] += reor;
        }    
    
        calculator_redfield calcR(evals);
        
        
    // setting up all transport rates and dephasings
    transport = 1;
    
    storage<asymptoticLF_complexv> g1un(1); // gfunctions first derivatives
    
    std::size_t found;
    found = experiment_complexity.find("ModRed");
    if(found!=std::string::npos)
	{
        if(!g1un.IsAlloc())
            g1un.Allocate(numBosc);
        for(int ind=0; ind<numBosc; ind++)
        {
            numericalSD spd(tempr,spdf.data1D[ind]);
            g1un.data1D[ind] = spd.GetGfD1f();
        }
        
        calcR.AddFluctuations(mfun,cfun,gfun,g1un,kij,gij,zijkl);
    }
    else
        calcR.AddFluctuations(mfun,kij,gij);
    

    
        // getting all reorganization energy matrix
        if(!reorganizations.IsAlloc())
        {
            reorganizations.Allocate(numT,numT);
            reorganizations = calcR.GetReorganizations();
        }


        
    storage<double> trates(2);
    if(found!=std::string::npos)
        trates = calcR.GetTransportRatesModRed();
    else
        trates = calcR.AddTransportRates();


    // GG block
    if(!ratesg.IsAlloc())
    {
        ratesg.Allocate(numG,numG);

        for(int ia=0; ia<numG; ia++)
        for(int ib=0; ib<numG; ib++)
            ratesg.data2D[ia][ib]=trates.data2D[ia][ib];
    }
            
    
    // EE block
    if(!ratese.IsAlloc())
    {
        ratese.Allocate(numE,numE);

        for(int ia=0; ia<numE; ia++)
        for(int ib=0; ib<numE; ib++)
            ratese.data2D[ia][ib]=trates.data2D[ia+numG][ib+numG];
    }    
    
    
        
        // compute dephasings
	if(!dephasings.IsAlloc())
	{
		if(nonmarkovian)
            dephasings = calcR.AddLifetimeDephasings();
        else
        {
            calcR.AddLifetimeDephasings();
            dephasings = calcR.AddPureDephasings();
        }
    }


	
	// set up the populations of the ground state
	if(!grPops.IsAlloc())
    {
        grPops.Allocate(numG);
        if(tempr == 0.0 || numG == 1)
        {
            if(numG>0)
                grPops.data1D[0] = 1.0;
        }
        else
        {
            // calculating the partition function
            double part = 0.0;
            double lowestE =evals.data1D[0]-reorganizations.data2D[0][0];
            for(int ind = 0; ind<numG; ind++)
                part += exp(-(evals.data1D[ind]-reorganizations.data2D[ind][ind]-lowestE)/tempr);
            for(int ind = 0; ind<numG; ind++)
                grPops.data1D[ind] = exp(-(evals.data1D[ind]-reorganizations.data2D[ind][ind]-lowestE)/tempr)/part;

            cout<<"# G band populations:\n";
            for(int ind = 0; ind<numG; ind++)
                cout<<grPops.data1D[ind]<<"\n";
        }
    }
}




