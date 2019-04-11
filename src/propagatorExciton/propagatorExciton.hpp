/* 
 * File:   propagatorExciton.hpp
 * Author: dariusa
 *
 * Created on October 29, 2015, 4:03 PM
 */

#ifndef PROPAGATOREXCITON_HPP
#define	PROPAGATOREXCITON_HPP

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


#include"../storage/storage.hpp"
#include"../complexv/complexv.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
//#include"../communicator_1st/communicator_1st.hpp"
#include"../communicator_3rd/communicator_3rd.hpp"
#include"../propagatorMemory/propagatorMemory.hpp"
#include"../calculator_redfield/calculator_redfield.hpp"

class propagatorExciton:
public propagatorMemory
{
public:
    
    // solves non-Markovian relaxation equation for Frenkel excitons
    // propagates 01, 11, 12 and 22 blocks
    
    // eigenstate basis
    
	propagatorExciton(int ig,int ie,int jf);
	//propagatorExciton(int ig); // only e block will be ready
        propagatorExciton(communicator_3rd& icom);
        //propagatorExciton(communicator_1st& icom);
    
	~propagatorExciton();
    // usage 
    // Initialize, Propagate and Exit
        
        
        storage<complexv> dmatrix0; // internal density matrix being propagated
        storage<complexv > dmatrixT; // density matrix with history
 
        void SetMarkovian(int);

        void SetEnergies(storage<double>& ien);
        void SetFluctuationMatrices(
		storage<double>& matr00,
		storage<double>& matr10,
		storage<double>& matr20,
		storage<double>& matr11,
		storage<double>& matr21,
		storage<double>& matr22);
        void SetFluctuationFunctions(storage<asymptoticLF_complexv>& imfun,storage<asymptoticLF_complexv>& icfun, double tempr);
    
        void MakeEquilibriumSim();
        
        // direct simulation routines
        void SetCurrentDM(storage<complexv>& imtr);
        
        //************************
        // nonmarkovian tools:

        // excitonic propagator when separately using cfun
        void ConvoluteCfun(storage<complexv>& der);
        void ConvoluteGen(storage<complexv>& der);
        storage<asymptoticLF_complexv> DMMcontinous;
        double DMMcontinousInternalTimeStep;
        int    DMMcontinousInternalTimeNump;
        storage<complexv> Propagate(storage<double>& times);



        // returns the whole history
        storage<complexv> GetHistory()
        {
            return dmatrixT;
        }
        void SetHistory(storage<complexv>& ihis)
        {
            // save the zero time initial condition
            storage<complexv> arr(1);
            int restoration = 0;
            int nt,ne;
            if(dmatrixT.IsAlloc())
            {
                dmatrixT.GetSize(nt,ne);
                arr.Allocate(ne);
                // saving zero time
                memcpy(arr.data1D,dmatrixT.data2D[0],ne*sizeof(complexv));
                restoration = 1;
            }
            dmatrixT = ihis;
            // restoring initial condition
            if(restoration)
                memcpy(dmatrixT.data2D[0],arr.data1D,ne*sizeof(complexv));
        }
        

        //***********************
        // generic functions
        storage<complexv>  PropagateDM(storage<double>& times);

        // sets the "delta" initial condition on a density matrix
	    void SetDeltaDM(int& ie,int &ig);

        // specifies the block of the excitonic density matrix
	    void SetBlock(int iblock);

public:
    
    int numG;
    int numE;
    int numF;
    
    int flagMarkovian;
    int flagMemoryWithCfun;
    int flagModred;
    int flagNonsecular;
    
    double meanEg;
    double meanEe;
    double meanEf;
    
    // eigenstate basis
    storage<double> energies;
 
    // resulting relaxation properties
    int superoperatorReady;
    storage<complexv> superoperatorSG; // secular markovian dephasings
    storage<double> superoperatorSP; // secular markovian populations
    storage<complexv> superoperatorM; // markovian full
    storage<complexv> superoperatorR; // non-markovian with memory
    storage<double> energiesReorg; // reorganization energies
    
    
    // redfield-type calculator for excitons
public:     calculator_redfield* calcR;
    
    
    
private:
    

	int block;
	int numL;
	int numR;
	int shL;
	int shR;
};

#endif	/* PROPAGATOREXCITON_HPP */

