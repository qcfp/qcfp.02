/* 
 * File:   propagatorExciton.hpp
 * Author: dariusa
 *
 * Created on October 29, 2015, 4:03 PM
 */

#ifndef PROPAGATOREXCITONWF_HPP
#define	PROPAGATOREXCITONWF_HPP

/* 
 * File:   propagatorMemory.hpp
 * Author: dariusa
 *
 * Created on October 29, 2015, 10:09 AM
 */

// this file defines
// class propagatorExciton (generic markovian and non-markovian (default) )
// class propagatorExcitonMarkovian (generic markovian)


#pragma once

#include"../toolsRandom/toolsRandom.hpp"
#include"../constants/constants.hpp"
#include"../storage/storage.hpp"
#include"../complexv/complexv.hpp"
#include"../communicator_3rd/communicator_3rd.hpp"
#include"../propagatorODEa/propagatorODEa.hpp"
class wavefunctionD1:
public communicator_3rd
{
public:
    
    storage<complexv> variationalparameters;
    double absolutetime;
    
    complexv* alpha0;
    // single-excitonic
    complexv* alpha1;
    // double-excitonic
    complexv** alpha2;
    
    
    int manifolds;
    int activemanifold;
    
    storage<int> manifold2sorter_u;
    storage<int> manifold2sorter_l;
    
    double meanEg;
    double meanEe;
    double meanEf;
    
    // a single spectral density is in use
    
    // bath
    complexv* lambda;
    
    int numQ;
    storage<double> bathfre;
    storage<double> bathsbc;
    
public:
    ~wavefunctionD1()
    {
        if(alpha2 != 0) delete[] alpha2;
    }
    
    wavefunctionD1()
    :communicator_3rd()
    {
        tempr = 0;
        variationalparameters.SetDimension(1);
        bathfre.SetDimension(1);
        bathsbc.SetDimension(2);
        absolutetime=0;
        manifold2sorter_u.SetDimension(1);
        manifold2sorter_l.SetDimension(1);
        alpha2 = 0;
        manifolds =0;
    }
    wavefunctionD1(int inb,int ine,int imanifolds)
    :communicator_3rd()
    {
        tempr = 0;
        variationalparameters.SetDimension(1);
        activemanifold = 0;
        numQ = inb;
        numE = ine;
        manifolds=imanifolds;
        manifold2sorter_u.SetDimension(1);
        manifold2sorter_l.SetDimension(1);
        absolutetime=0;
        alpha2 = 0;
        if(manifolds==1)
        {
            variationalparameters.Allocate(numQ+numE);
            alpha0 = 0;
        }
        if(manifolds==2)
        {
            variationalparameters.Allocate(numQ+numE+1);
            alpha0 = variationalparameters.data1D+numQ+numE;
        }
        else if(manifolds==3)
        {
            variationalparameters.Allocate(numQ+numE+1+numE*numE);
            alpha0 = variationalparameters.data1D+numQ+numE;
            alpha2 = new complexv*[numE];
            for(int ie=0;ie<numE;ie++)
                alpha2[ie] = variationalparameters.data1D+numQ+numE+1 + numE*ie;
        }
        lambda = variationalparameters.data1D;
        alpha1 = variationalparameters.data1D+numQ;
        bathfre.SetDimension(1);
        bathsbc.SetDimension(2);
        
        bathfre.Allocate(numQ);
        bathsbc.Allocate(numE,numQ);
    }
    
    
    
    
    
    
    
    wavefunctionD1(communicator_3rd& comm)
    :communicator_3rd()
    {
        /*
         *		if(!comm.excitonictransformed)
         *		{
         *			std::cout<<"Error: please transform to eigenstates before calling D wavefunctions\n";
         *			return;
    }
    excitonictransformed=comm.excitonictransformed;
    
    
    if(comm.coupl.CheckDimension()==3)
    {
    std::cout<<"Error: this kind of correlations is not supported; instead use many spectral density functions\n";
    return;
    }
    */
        
        inputfile<<comm.inputfile.str();
        tempr = comm.tempr;
        
        absolutetime=0;
        numG = comm.numG;
        numE = comm.numE;
        numF = comm.numF;
        
        band = comm.band; 
        // 1 - one exciton properties
        // 2 - 1 and 2 exciton properties 
        
        // constructing bath:
        if(true){
            // first reading the vibrational modes from the input file
            toolsIO tio;
            string str = "";
            storage<double> arr(2);
            storage<int> iri(2);
            int nummodes = 0;
            // setting parameters from input communicator_3rd
            // list of specific vibrational frequencies and their couplings
            if(tio.LookUpAndReadSquareMatrix<int>(
                "BathSpecificModes:",
                "Reading frequencies and coupling strengths of specific bath modes. Coupling term is |n><n|*w*c*q. The reorganization energy of that mode is L=w(c^2)/2, HR parameter is S=(c^2)/2. \n", 
                                                     1,1, iri,inputfile.str()))
            {
                // first reading the number of modes and then the whole list of modes
                nummodes= iri.data2D[0][0];
                
                // now reading the list of frequencies and their couplings to distinct sites
                string sstr = inputfile.str();
                string str = "";
                // looking for entry
                std::size_t found = sstr.find("BathSpecificModes:");
                sstr = sstr.substr(found);
                std::istringstream ifs(sstr);
                getline(ifs,str);
                getline(ifs,str);
                if(arr.IsAlloc()) arr.Delete();
                arr.Allocate(nummodes,numE+1);
                if(tio.ReadRectangular(arr, ifs))
                {
                    cout<<"Error: Some problem with file reading\n";
                }		// assignment
            }
            // next counting the spectral density discretized baths
            // setting the specific discretized bath
            int numtot=0;
            int numbo=comm.spdf.GetSize();
            for(int ib=0;ib<numbo;ib++)
            {
                numtot+=comm.spdf.data1D[ib].GetN();
            }
            // so the total number of modes becomes equal to
            numQ=numtot+nummodes;
            
            bathfre.SetDimension(1);
            bathsbc.SetDimension(2);
            bathfre.Allocate(numQ);
            bathsbc.Allocate(numE,numQ);
            
            
            for(int ind=0; ind<nummodes; ind ++){
                bathfre.data1D[ind]=arr.data2D[ind][0];
                for(int ins=0; ins<numE; ins ++){
                    bathsbc.data2D[ins][ind]=arr.data2D[ind][ins+1];
                }
            }
            numtot = 0;
            for(int ib=0;ib<numbo;ib++)
            {
                int nump = comm.spdf.data1D[ib].GetN();
                double fref=comm.spdf.data1D[ib].GetXF();
                double fres=comm.spdf.data1D[ib].GetStep();
                double* value=comm.spdf.data1D[ib].DirectAccessD();
                
                for(int inv=0;inv<nump;inv++)
                {
                    double fff = fres*inv;
                    bathfre.data1D[numtot+nummodes]=fff;
                    
                    for(int ins=0;ins<numE;ins++)
                        if(inv!=0)
                            bathsbc.data2D[ins][numtot+nummodes]=sqrt( (fres/fff)*(comm.coupl.data2D[ib][ins])*(value[inv]/fff)/constants::piover2);
                        numtot +=1;
                }
            }
        }
        
        
        // 
        manifold2sorter_u.SetDimension(1);
        manifold2sorter_l.SetDimension(1);
        
        manifolds = 0;
        if(band == 1)
        {
            manifolds = 2;
            activemanifold = 0;
        }
        else if(band == 2)
        {
            manifolds = 3;
            activemanifold = 0;
        }
        bosonic = comm.bosonic;
        
        // the comm object must contain all excitonic information
        
        
        /////////
        // excitonic system:
        bosonic = comm.bosonic;
        ham=comm.ham;
        anh = comm.anh;
        dip = comm.dip;
        dip2Corrections = comm.dip2Corrections;
        coupl = comm.coupl;
        
        edips = comm.edips;
        
        evec0 = comm.evec0;
        evec1 = comm.evec1;
        evec2 = comm.evec2;
        
        // shifting hamiltonian by the vibrational reorganization energies:
        for(int id=0;id<numE;id++)
        {
            double reoe = 0.0;
            for(int iq=0;iq<numQ;iq++)
                reoe += bathfre.data1D[iq]*bathsbc.data2D[id][iq]*bathsbc.data2D[id][iq]/2.0;

            std::cout<<"Shifting hamiltonian and anharmonicities by reorganization energies:\n";
            std::cout<<id<<" "<<reoe<<"\n";
            ham.data2D[id][id] += reoe;
            
            if(anh.IsAlloc())
                anh.data2D[id][id] += 4*reoe;
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
        // mean energies
        
        meanEg = 0.0;
        meanEe = 0.0;
        meanEf = 0.0;
        for(int id=0;id<numE;id++)
        {
            meanEe += ham.data2D[id][id];
            
            if(band == 2){
                if(bosonic)
                {
                    for(int ig=0;ig<numE;ig++)
                        meanEf += ham.data2D[id][id]+ham.data2D[ig][ig]+anh.data2D[id][ig];
                }
                else
                {
                    for(int ig=0;ig<numE;ig++)
                        if(ig!=id)
                            meanEf += ham.data2D[id][id]+ham.data2D[ig][ig]+anh.data2D[id][ig];
                }
            }
        }
        meanEe /= numE;
        meanEf /= numE*(numE-1);
        
        
        variationalparameters.SetDimension(1);
        
        alpha2 = 0;
        if(manifolds==1)
        {
            variationalparameters.Allocate(numQ+numE);
            alpha0 = 0;
        }
        if(manifolds==2)
        {
            variationalparameters.Allocate(numQ+numE+1);
            alpha0 = variationalparameters.data1D+numQ+numE;
        }
        else if(manifolds==3)
        {
            variationalparameters.Allocate(numQ+numE+1+numE*numE);
            alpha0 = variationalparameters.data1D+numQ+numE;
            alpha2 = new complexv*[numE];
            for(int ie=0;ie<numE;ie++)
                alpha2[ie] = variationalparameters.data1D+numQ+numE+1 + numE*ie;
            
            //			manifold2sorter_u.Allocate(numF);
            //			manifold2sorter_l.Allocate(numF);
            //			int tind=0;
            //			for(int iu=0;iu<numE;iu++)
            //			{
            //				int numl = (bosonic)?iu+1:iu;
            //				for(int il=0;il<numl;il++)
            //				{
            //					manifold2sorter_u.data1D[tind]=iu;
            //					manifold2sorter_l.data1D[tind]=il;
            //					tind++;
            //				}
            //			}
            
            manifold2sorter_u = comm.ham2sorter_l;
            manifold2sorter_l = comm.ham2sorter_r;
        }
        lambda = variationalparameters.data1D;
        alpha1 = variationalparameters.data1D+numQ;
        
    }
    
    
    /////////////////////////////////////////////////////////
    void SetLambda(storage<complexv> ilam)
    {
        if(numQ != ilam.GetSize())
        {
            cout<<"Error in wavefunctionD1:SetVibrations: sizes do not match\n";
            return;
        }
        if(!variationalparameters.IsAlloc())
        {
            cout<<"Error in wavefunctionD1:SetVibrations: destination is not allocated\n";
            return;
        }
        memcpy(lambda,ilam.data1D,numQ*sizeof(complexv));
        absolutetime=0;
    }
    void SetAlpha1(storage<complexv> ialp)
    {
        if(numE != ialp.GetSize())
        {
            cout<<"Error in wavefunctionD1:SetVibrations: sizes do not match\n";
            return;
        }
        if(!variationalparameters.IsAlloc())
        {
            cout<<"Error in wavefunctionD1:SetVibrations: destination is not allocated\n";
            return;
        }
        memcpy(alpha1,ialp.data1D,numE*sizeof(complexv));
        absolutetime=0;
    }
    void CopyWF(wavefunctionD1& iwf)
    {
        variationalparameters.Delete();
        if(alpha2!=0) delete[] alpha2;
        bathfre.Delete();
        bathsbc.Delete();
        variationalparameters = iwf.variationalparameters;
        activemanifold = iwf.activemanifold;
        manifolds = iwf.manifolds;
        numQ = iwf.numQ;
        numE = iwf.numE;
        numF = iwf.numF;
        evec0 = iwf.evec0;
        evec1 = iwf.evec1;
        evec2 = iwf.evec2;
        edips = iwf.edips;
        
        ham = iwf.ham;
        anh = iwf.anh;
        dip = iwf.dip;
        dip2Corrections = iwf.dip2Corrections;
        
        meanEg = iwf.meanEg;
        meanEe = iwf.meanEe;
        meanEf = iwf.meanEf;
        
        bathfre=iwf.bathfre;
        bathsbc=iwf.bathsbc;
        
        alpha2 = 0;
        if(manifolds==1)
        {
            alpha0 = 0;
        }
        if(manifolds==2)
        {
            alpha0 = variationalparameters.data1D+numQ+numE;
        }
        else if(manifolds==3)
        {
            alpha0 = variationalparameters.data1D+numQ+numE;
            alpha2 = new complexv*[numE];
            for(int ie=0;ie<numE;ie++)
                alpha2[ie] = variationalparameters.data1D+numQ+numE+1 + numE*ie;
        }
        lambda = variationalparameters.data1D;
        alpha1 = variationalparameters.data1D+numQ;
        
        absolutetime=iwf.absolutetime;
    }
    
    void OpticalTransition01(int ini,int ifi)
    {
        activemanifold = 1;
        //result =0;
        // in the present model there is only one electronic ground state
        // optical transition happens in exciton eigenbasis
        //for(int ina=0;ina<numE;ina++)
        //	alpha1[ina]=evec1.data2D[ifi][ina]*alpha0[0];
        
        memset(alpha1,0,numE*sizeof(complexv));
        alpha1[ifi]=alpha0[0];
        //return dip.data1D[ifi];
        
    }
    void OpticalTransition10(int ini,int ifi)
    {
        activemanifold = 0;
        //result =0;
        // in the present model there is only one electronic ground state
        // optical transition happens in exciton eigenbasis
        //alpha0[0]=0;
        //for(int ina=0;ina<numE;ina++)
        //	alpha0[0]+=evec1.data2D[ini][ina]*alpha1[ina];
        
        alpha0[0]=alpha1[ini];
        //return dip.data1D[ini];
    }
    
    void OpticalTransition12(int ini,int ifi)
    {
        activemanifold = 2;
        //result =0;
        // in the present model there is only one electronic ground state
        // optical transition happens in exciton eigenbasis
        //for(int ina=0;ina<numE;ina++)
        //for(int inb=0;inb<numE;inb++)
        //	alpha2[ina][inb]=evec2.data3D[ifi][ina][inb]*evec1.data2D[ini][ina]*alpha1[ina];
        
        memset(alpha2[0],0,numE*numE*sizeof(complexv));
        
        int ifi1=manifold2sorter_u.data1D[ifi];
        int ifi2=manifold2sorter_l.data1D[ifi];
        
        if(ifi1 == ini && ifi2 != ini)
        {
            alpha2[ifi1][ifi2]=alpha1[ini];
            alpha2[ifi2][ifi1]=alpha1[ini];
            return;// dip.data1D[ifi2];
        }
        if(ifi2 == ini && ifi1 != ini)
        {
            alpha2[ifi1][ifi2]=alpha1[ini];
            alpha2[ifi2][ifi1]=alpha1[ini];
            return;// dip.data1D[ifi1];
        }
        
        if(bosonic)
        {
            if(ifi2 == ini && ifi1 == ini)
            {
                alpha2[ifi1][ifi2]=alpha1[ini];
                return;// (dip.data1D[ini]+ ((dip2Corrections.IsAlloc())?(dip2Corrections.data1D[ini]):dvector3d(0,0,0)) )*constants::root2;
            }
            else
            {
                //result = 1;
                return ;//dvector3d(0,0,0);
            }
        }
        else
        {
            //result = 1;
            return;// dvector3d(0,0,0);
        }
    }
    void OpticalTransition21(int ini,int ifi)
    {
        activemanifold = 1;
        //result =0;
        memset(alpha1,0,numE*sizeof(complexv));
        
        int ini1=manifold2sorter_u.data1D[ini];
        int ini2=manifold2sorter_l.data1D[ini];
        
        if(ini1 == ifi && ini2 != ifi)
        {
            alpha1[ifi]=alpha2[ini1][ini2];
            return;// dip.data1D[ini2];
        }
        if(ini2 == ifi && ini1 != ifi)
        {
            alpha1[ifi]=alpha2[ini1][ini2];
            return;// dip.data1D[ini1];
        }
        if(bosonic)
        {
            if(ini2 == ifi && ini1 == ifi)
            {
                alpha1[ifi]=alpha2[ini1][ini2];
                return ;//(dip.data1D[ifi]+ ((dip2Corrections.IsAlloc())?dip2Corrections.data1D[ifi]:dvector3d(0,0,0)) )*constants::root2;
            }
            else
            {
                //result = 1;
                return ;//dvector3d(0,0,0);
            }
        }
        else
        {
            //result = 1;
            return ;//dvector3d(0,0,0);
        }
    }
    
    // from frequencies creates initial shifts
    void MakeThermal()
    {
        absolutetime = 0;
        toolsRandom rnd;
        //        int numb = cfun.GetN();
        //        double fini = cfun.GetXI();
        //        double fste = cfun.GetStep();
        
        activemanifold = 0;
        for(int i0=0;i0<numE;i0++)
            alpha1[i0]=0.0;
        
        if(manifolds==2)
            *alpha0 = 1.0;
        
        if(manifolds == 3)
        {
            *alpha0 = 1.0;
            for(int i0=0;i0<numE;i0++)
                for(int i1=0;i1<numE;i1++)
                    alpha2[i0][i1]=0.0;
        }
        
        // using Glauber Sudarshan representation
        //std::cout<<"Error: thermal WF needs to be rewritten\n";
        for(int i0=0;i0<numQ;i0++)
        {
            double freq = bathfre.data1D[i0];
            
            if(tempr == 0.0)
            {
                lambda[i0] = 0.0;
            }
            else
            {
                if(freq>constants::smallepsilon)
                {
                    double broad = 0.5/sqrt(exp(freq/tempr)-1);
                    lambda[i0] = rnd.RandGD(0,broad);
                    //lambda[i0] = rnd.RandED(sqrt(2.0*tempr/freq/freq));
                }
                else
                    lambda[i0] = 0.0;
            }
            lambda[i0] *= exp(coni*rnd.RandLD(constants::pi2));
        }
    }
    
    
    void CopyState(wavefunctionD1& iwf)
    {
        absolutetime = iwf.absolutetime;
        if(numQ != iwf.numQ || numE != iwf.numE)
        {
            cout<<"Error in CopyState0(wavefunctionD1 iwf)\n";
            return;
        }
        memcpy(lambda,iwf.lambda,numQ*sizeof(complexv));
        if(iwf.activemanifold == 0)
        {
            *alpha0 = *(iwf.alpha0);
            activemanifold =0;
        }
        else if(iwf.activemanifold == 1)
        {	
            memcpy(alpha1,iwf.alpha1,numE*sizeof(complexv));
            activemanifold = 1;
        }
        else if(iwf.activemanifold == 2)
        {
            std::cout<<"Error in WF: cannot copy WF of the 2-nd manifold: code is not ready\n";
        }
        return;
    }
    
    complexv GetBathOverlap(wavefunctionD1& iwf)
    {
        complexv retval = 0.0;
        for(int ind=0;ind<numQ;ind++)
            retval += -0.5* (lambda[ind].norm()+iwf.lambda[ind].norm()-2.0*lambda[ind]*iwf.lambda[ind].conjugate()* exp(cnni*bathfre.data1D[ind]*(absolutetime-iwf.absolutetime)));
        
        return exp(retval);
    }
    complexv GetQuantumOverlap(wavefunctionD1& iwf)
    {
        complexv retval = 0.0;
        
        if(activemanifold!=iwf.activemanifold)
            return 0;
        if(activemanifold == 0)
            return alpha0[0]*(iwf.alpha0[0]).conjugate();
        
        if(activemanifold == 1)
            for(int ind=0; ind<numE;ind++)
                retval += alpha1[ind]*(iwf.alpha1[ind]).conjugate();
            
            if(activemanifold == 2)
                for(int ind1=0; ind1<numE;ind1++)
                    for(int ind2=0; ind2<=ind1;ind2++)
                        retval += alpha2[ind1][ind2]*(iwf.alpha2[ind1][ind2]).conjugate();
                    
                    
                    return retval;
    }
    
    void operator ()(double itime, complexv* state_now,  complexv* derivatives)
    // returns an array of derivative values based on the ODE
    {
        
        if(activemanifold == 1)
        {
            complexv* alpha=state_now+numQ;
            complexv* lambda=state_now;
            complexv* dalpha=derivatives+numQ;
            complexv* dlambda=derivatives;
            
            
            double** ham2D = ham.data2D;
            double** gval2D = bathsbc.data2D;
            double* bfre1D = bathfre.data1D;
            for(int ind=0;ind<numE;ind++)
            {
                dalpha[ind]= cnni*alpha[ind]*(ham2D[ind][ind]-meanEe);
                for(int inp=0;inp<ind;inp++)
                    dalpha[ind] += cnni*alpha[inp]*ham2D[ind][inp];
                for(int inp=ind+1;inp<numE;inp++)
                    dalpha[ind] += cnni*alpha[inp]*ham2D[ind][inp];
                
                
                
                for(int inq=0;inq<numQ;inq++)
                    dalpha[ind] += coni*alpha[ind]*gval2D[ind][inq]/constants::root2*bfre1D[inq]*(lambda[inq]*exp(cnni*bfre1D[inq]*itime)).real();
//                    dalpha[ind] += coni*gval2D[ind][inq]*(lambda[inq]*exp(cnni*bfre1D[inq]*itime)).real();
                
//                cout<<itime<<" "<<dalpha[0].re()<<" "<<dalpha[0].im()<<" "<<(lambda[0]*exp(cnni*bfre1D[0]*itime)).real()<<" "<<(lambda[0]*exp(cnni*bfre1D[0]*itime)).imag()<<"\n";
            }
            for(int inq=0;inq<numQ;inq++)
            {
                dlambda[inq]=0.0;//cnni*bfre1D[inq]*lambda[inq];
                for(int inp=0;inp<numE;inp++)
                    dlambda[inq]+= exp(coni*bfre1D[inq]*itime)* coni*bfre1D[inq]*gval2D[inp][inq]/constants::root2*(alpha[inp]).norm();
            }
            
        }
        else if(activemanifold == 2)
        {
            complexv* alpha=state_now+numQ;
            complexv* lambda=state_now;
            complexv* dalpha=derivatives+numQ;
            complexv* dlambda=derivatives;
            
            
            double** ham2D = ham.data2D;
            double** gval2D = bathsbc.data2D;
            double* bfre1D = bathfre.data1D;
            
            int addanh = (anh.IsAlloc())? 1: 0;
            
            for(int ind=0;ind<numE;ind++){
                int addupto = (bosonic)? ind+1 : ind; 
                
                for(int inr=0;inr<addupto;inr++){
                    
                    // diagonal part
                    dalpha[ind*numE+inr]= cnni*alpha[ind*numE+inr]*(ham2D[ind][ind]+ham2D[inr][inr]-meanEf);
                    
                    if(addanh)
                        dalpha[ind*numE+inr]= cnni*alpha[ind*numE+inr]*anh.data2D[ind][inr];
                    
                    // simple couplings
                    if(ind != inr){
                        // here I have combination band
                        for(int inp=0;inp<inr;inp++)
                            dalpha[ind*numE+inr] += cnni*alpha[ind*numE+inp]*ham2D[inr][inp];
                        for(int inp=inr+1;inp<numE;inp++)
                            dalpha[ind*numE+inr] += cnni*alpha[ind*numE+inp]*ham2D[inr][inp];
                        
                        for(int inp=0;inp<ind;inp++)
                            dalpha[ind*numE+inr] += cnni*alpha[inp*numE+inr]*ham2D[ind][inp];
                        for(int inp=ind+1;inp<numE;inp++)
                            dalpha[ind*numE+inr] += cnni*alpha[inp*numE+inr]*ham2D[ind][inp];
                        
                        // overtone coupling corrections
                        dalpha[ind*numE+inr] += cnni*alpha[ind*numE+ind]*(ham2D[inr][ind]*(constants::root2-1));
                        dalpha[ind*numE+inr] += cnni*alpha[inr*numE+inr]*(ham2D[ind][inr]*(constants::root2-1));
                    }else{
                        // here ind == inr -- overtones
                        for(int inp=0;inp<ind;inp++)
                            dalpha[ind*numE+ind] += cnni*alpha[ind*numE+inp]*ham2D[ind][inp]*constants::root2;
                        for(int inp=ind+1;inp<numE;inp++)
                            dalpha[ind*numE+ind] += cnni*alpha[ind*numE+inp]*ham2D[ind][inp]*constants::root2;
                        
                    }
                    
                    // couplings to the bath
                    for(int inq=0;inq<numQ;inq++)
                        dalpha[ind+numE*inr] += coni*alpha[ind*numE+inr]*(gval2D[ind][inq]+gval2D[inr][inq])/constants::root2*bfre1D[inq]*(lambda[inq]*exp(cnni*bfre1D[inq]*itime)).real();
                    
                    dalpha[inr+numE*ind]=dalpha[ind+numE*inr];
                    
                }}
                for(int inq=0;inq<numQ;inq++)
                {
                    dlambda[inq]=0.0;//cnni*bfre1D[inq]*lambda[inq];
                    for(int inp=0;inp<numE;inp++)
                        for(int inr=0;inr<=inp;inr++)
                            dlambda[inq]+= exp(coni*bfre1D[inq]*itime)* coni*bfre1D[inq]*(gval2D[inp][inq]+gval2D[inr][inq])/constants::root2*(alpha[inp*numE+inr]).norm();
                }
        }
    }
    
};


class propagatorExcitonWF
{
    // solves non-Markovian relaxation equation for Frenkel excitons
    // propagates 01, 11, 12 and 22 blocks
    
public:
    double tstep;
    
    private: propagatorODEa<complexv,wavefunctionD1>* objprop;
    
public:
    
    
    propagatorExcitonWF(double istep)
    {
        tstep = istep;
        objprop = 0;
    }
    
    ~propagatorExcitonWF()
    {
        if( objprop != 0) delete objprop;
        objprop = 0;
    }
    // usage 
    // Initialize, Propagate and Exit
    
    void TransferState01(int ie,wavefunctionD1& iwf)
    {
        iwf.activemanifold = 1;
        iwf.alpha1[ie] = *iwf.alpha0;
    }
    void TransferState10(int ie,wavefunctionD1& iwf)
    {
        iwf.activemanifold = 0;
        *iwf.alpha0 = iwf.alpha1[ie];
    }
    void TransferState12(int ie,int iz,wavefunctionD1& iwf)
    {// from ie to ie_iz
        iwf.activemanifold = 2;
        iwf.alpha2[ie][iz] = iwf.alpha1[ie];
        iwf.alpha2[iz][ie] = iwf.alpha1[ie];
    }
    void TransferState21(int iz,int ie,wavefunctionD1& iwf)
    {// from iz_ie to ie
        iwf.activemanifold = 1;
        iwf.alpha1[ie] = iwf.alpha2[ie][iz];
    }
    
    void SwitchManifold(int imani,wavefunctionD1& iwf)
    {
        if(iwf.manifolds == 1)
        {
            if(imani !=1)
            {
                std::cout<<"Error: manifold is incorrect in SwitchManifold()\n";
                return;
            }
        }
        else if(imani >= iwf.manifolds)
        {
            std::cout<<"Error: manifold is incorrect in SwitchManifold()\n";
            return;
        }
        iwf.activemanifold = imani;
    }
    
    void PropagateWF(double time,wavefunctionD1& iwf)
    {
        if(iwf.activemanifold == 0)
        {
            // propagate ground state
            // electronic amplitude stays constant
            
            // nuclear parameter performs rotation around zero
            // it is constant in rwa approach
            //for(int ib=0; ib<iwf.numQ; ib++)
            //{
            //    iwf.lambda[ib] = exp(cnni*iwf.bathfre.data1D[ib]*(time-iwf.absolutetime))*iwf.lambda[ib];
           // }
            iwf.absolutetime=time;
        }
        else if(iwf.activemanifold > 0 )
        {
            // propagate single excited state
            // electronic amplitude rotates
            
            // propagation with small internal timesteps;
            int numloops = (int)((time-iwf.absolutetime)/tstep); // this should return floor value
            double reminder = time-iwf.absolutetime-numloops*tstep;
            
            if(reminder<constants::smallepsilon)
            {
                if(numloops>0)
                {
                    numloops--;
                    reminder += tstep;
                }
                else
                {
                    //skip calculation
                    return;
                }
            }
            
            // number of degrees of freedom
            int numbtotal;
            if(iwf.activemanifold == 1)
                numbtotal = iwf.numE+iwf.numQ;
            else
                numbtotal = iwf.numE*iwf.numE+iwf.numQ;
            
            if(objprop == 0) objprop = new propagatorODEa<complexv,wavefunctionD1>(&iwf, numbtotal);
            else if(objprop->GetNumDF() != numbtotal)
            {
                delete objprop;
                objprop = new propagatorODEa<complexv,wavefunctionD1>(&iwf, numbtotal);
            }
            
            storage<complexv> state(1);
            state.Allocate(numbtotal);
            complexv* state_current=state.data1D;
            
            if(iwf.activemanifold == 1)
                memcpy(state_current,iwf.variationalparameters.data1D,numbtotal*sizeof(complexv));
            else
            {
                memcpy(state_current, iwf.variationalparameters.data1D, iwf.numQ*sizeof(complexv));
                memcpy(state_current+iwf.numQ, iwf.variationalparameters.data1D+iwf.numQ+iwf.numE+1, iwf.numE*iwf.numE*sizeof(complexv));
            }
            
            
            for(int it=0;it<numloops;it++)
            {
                int solution = objprop->propagateODEstep(iwf.absolutetime+it*tstep, state_current, tstep);
            }
            int solution = objprop->propagateODEstep(iwf.absolutetime+numloops*tstep, state_current, reminder);
            
            // updating wavefunction
            iwf.absolutetime = time;
            
            if(iwf.activemanifold == 1)
                memcpy(iwf.variationalparameters.data1D,state_current,numbtotal*sizeof(complexv));
            else
            {
                memcpy(iwf.variationalparameters.data1D, state_current, iwf.numQ*sizeof(complexv));
                memcpy(iwf.variationalparameters.data1D+iwf.numQ+iwf.numE+1, state_current+iwf.numQ, iwf.numE*iwf.numE*sizeof(complexv));
            }
        }
        
        
        return;
    }
    
    
    
    
};


#endif	/* PROPAGATOREXCITON_HPP */

