// This project provides means for coomputing orientationally
// averaged transition amplitudes.
// Two and four vectors are included
// double vector dvector3d class is used for fast operations

#pragma once
#include"../dvector3d/dvector3d.hpp"
#include"../dtensor3x3/dtensor3x3.hpp"
#include"../complexv/complexv.hpp"


class averages_orient
{
public:
    
	// this constructor includes setting the type
	// of averaging:
	// 0 - isotropic averaging
	// 1 - no averaging (explicit products):
    // in this case local and lab frames coincide:
    averages_orient(int itype);
    averages_orient();
    
    // additionally sets the type of averaging
	// 0 - isotropic averaging
	// 1 - no averaging (explicit products):
    // in this case local and lab frames coincide:
    void SetType(int itype);

    // a member function provides the averaging of two vectors:
    // es[0] and es[1] are field vectors in the lab frame
    // ds[0] and ds[1] are system vectors in the local frame
    double rot_av_dip_2(dvector3d* es,dvector3d* ds);

    // a member function provides the averaging of four vectors:
    // es[0] ... to es[3] are field vectors in the lab frame
    // ds[0] ... to ds[3] are system vectors in the local frame
	double rot_av_dip_4(dvector3d* es,dvector3d* ds);



	complexv rot_av_complete_2_4(
		 dvector3d* es, // optical polarizations
		 dvector3d* ks, // optical wavevectors
		 double* om, // fourier frequencies
		 dvector3d* ds, // electric transition dipoles
		 dvector3d* ms, // magnetic transition dipoles
		 dtensor3x3* ts, // quadrupole transition tensors
		 int number // number of interactions
		 );

private:
    
	// the local type variable
    int type;

    // a member function provides the isotropic averaging of two vectors:
    // es[0] and es[1] are field vectors in the lab frame
    // ds[0] and ds[1] are system vectors in the local frame
    double rot_av_dip_iso_2(dvector3d* es,dvector3d* ds);

    // a member function provides the explicit amplitude of two vectors
    // local and lab frames coincide:
    // es[0] and es[1] are field vectors in the lab frame
    // ds[0] and ds[1] are system vectors in the local frame
    double rot_av_dip_exp_2(dvector3d* es,dvector3d* ds);

    // a member function provides the isotropic averaging of four vectors:
    // es[0] ... to es[3] are field vectors in the lab frame
    // ds[0] ... to ds[3] are system vectors in the local frame
    double rot_av_dip_iso_4(dvector3d* es,dvector3d* ds);

    // a member function provides the explicite amplitude of four vectors:
    // local and lab frames coincide:
    // es[0] ... to es[3] are field vectors in the lab frame
    // ds[0] ... to ds[3] are system vectors in the local frame
    double rot_av_dip_exp_4(dvector3d* es,dvector3d* ds);


};

