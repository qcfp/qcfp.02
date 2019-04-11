#pragma once
#include"../complexv/complexv.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
#include"../propagatorM/propagatorM.hpp"
#include"../correlation4points/correlation4points.hpp"


// diagram
        
//       e |  | e
//        -------
//       k |  | l
//       d |  | f
//        -------
//       j |  | m
//       c |  | g
//        -------
//       i |  | n
//       b |  | h
//        -------
//       a |  | a

//                states[4 ] /*|*/ states[4 ] 
//                //-------------------------
//                states[10] /*|*/ states[11] 
//                states[3 ] /*|*/ states[5 ] 
//                //-------------------------
//                states[9 ] /*|*/ states[12] 
//                states[2 ] /*|*/ states[6 ] 
//                //-------------------------
//                states[8 ] /*|*/ states[13] 
//                states[1 ] /*|*/ states[7 ] 
//                //-------------------------
//                states[0 ] /*|*/ states[0 ] 



// input techniques:
//               technique ==1 // llll
//               technique ==2 // lllr
//               technique ==3 // llrl
//               technique ==4 // llrr
//               technique ==5 // lrll
//               technique ==6 // lrlr
//               technique ==7 // lrrl
//               technique ==8 // lrrrr
//    
//        // incoherent diagrams
//               technique ==9 // kii 
//               technique ==10 // ki


class feinman2sideddiagram3
{
private:

    
        int* states;//[14];


    // energy (imaginary part is for the dephasing)
	complexv* freq; //freq[3];

	// amplitudes of transitions (dipoles)
	//dvector3d *ds;

	// vectors of fields (dipoles)
	//dvector3d *es;

        // amplitude of the diagram
        complexv damplitude;

        // bath:
        
	asymptoticLF<double>* grF;
	propagatorM* propF;
	int ongrF;
        
//        double*** gijAmplitude;
//        asymptoticLF_complexv** gij;
	asymptoticLF_complexv *gf;
	double ***gfamp;
        int technique;
        correlation4points* c4p_obj;
        int inumO;
	int ongij;
        double ***data_real, ***data_imag, *freq_real, *freq_imag, *gf_real, *gf_imag;
        asymptoticLF_double *gRe;
	asymptoticLF_double *gIm;
        storage<double>* gRdat, *grdat;
	storage<double>* gIdat;
        double *gRD, *grd;
	double *gID;
        int nstates, taud, tauc, taub, taua;
				    
        int flagCheck;

        
public:
    
    double feinman2sideddiagram3smallparameter;
    
feinman2sideddiagram3();
feinman2sideddiagram3(complexv* ifreq, dvector3d* ids, dvector3d* ies,int& averaging);
feinman2sideddiagram3(complexv* ifreq, complexv factor);
~feinman2sideddiagram3();

void assignLineshape(int& itechnique,asymptoticLF_complexv* igij,double*** igijAmplitude,int& inumOsc,int* istates, int&);
void assignPropagator(asymptoticLF<double>* iGij);
void assignPropagator(propagatorM* iGij);
//void multiplier(double&);
void propagate(storage<complexv>& res,double* itT, int* itN);
complexv calculate(double&,double&,double&);
//void printout();

};

