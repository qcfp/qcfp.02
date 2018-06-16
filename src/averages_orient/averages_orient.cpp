#include<iostream>
#include<cstring>
#include"averages_orient.hpp"
#include"../constants/constants.hpp"

averages_orient::averages_orient()
{
    type = 0;
}

averages_orient::averages_orient(int itype)
{
    type = itype;
}

void averages_orient::SetType(int itype)
{
    type = itype;
}

double averages_orient::rot_av_dip_2(dvector3d* es,dvector3d* ds)
{
    if(type == 0)
        return rot_av_dip_iso_2(es,ds);
    else
        return rot_av_dip_exp_2(es,ds);
}


double averages_orient::rot_av_dip_4(dvector3d* es,dvector3d* ds)
{
    if(type == 0)
        return rot_av_dip_iso_4(es,ds);
    else
        return rot_av_dip_exp_4(es,ds);
}


double averages_orient::rot_av_dip_iso_2(dvector3d* es,dvector3d* ds)
{
    
        // two vectors averaged
	// DIPOLE APPROXIMATION
			
	double averaging = (es[1]*es[0]) * (ds[1]*ds[0]);
	return averaging / 3.0;
}

double averages_orient::rot_av_dip_exp_2(dvector3d* es,dvector3d* ds)
{
        // two vectors explicitly used
	// DIPOLE APPROXIMATION
			
	double averaging = (es[1]*ds[1]) * (es[0]*ds[0]);
	return averaging;
}

double averages_orient::rot_av_dip_iso_4(dvector3d* es,dvector3d* ds)
{


    double averaging=0;

    int mm[3][3] = {
		{ 4, -1, -1},
		{-1,  4, -1},
		{-1, -1,  4}
	};

	double rivec[3];
	double levec[3];
	double revec[3];
    
        
        // Getting values:
	rivec[0] = (ds[3]*ds[2])*(ds[1]*ds[0]);
	rivec[1] = (ds[3]*ds[1])*(ds[2]*ds[0]);
	rivec[2] = (ds[3]*ds[0])*(ds[2]*ds[1]);
			
	levec[0]  = (es[3]*es[2])*(es[1]*es[0]);
	levec[1]  = (es[3]*es[1])*(es[2]*es[0]);
	levec[2]  = (es[3]*es[0])*(es[2]*es[1]);
			
	// 1 right vector times middle matrix
	for(int indi=0;indi<3;indi++)
	{
		revec[indi]=0;
		for(int indk=0;indk<3;indk++)
			revec[indi] += mm[indi][indk]*rivec[indk];
	}
            
	// 2 lvector times result_vector
	averaging=0;
	for(int indk=0;indk<3;indk++)
		averaging += levec[indk]*revec[indk];
			
	averaging /= 30.0;
        return averaging;
      
}

double averages_orient::rot_av_dip_exp_4(dvector3d* es,dvector3d* ds)
{
        // two vectors explicitly used
	// DIPOLE APPROXIMATION
			
	double averaging = (es[3]*ds[3]) * (es[2]*ds[2]) * (es[1]*ds[1]) * (es[0]*ds[0]);
	return averaging;
}





// orientational average for an arbitrary response up to third order
complexv averages_orient::rot_av_complete_2_4(
				 dvector3d* es, // optical polarizations
				 dvector3d* ks, // optical wavevectors
				double* om, // entering Fourier frequencies
				 dvector3d* ds, // electric transition dipoles
				 dvector3d* ms, // magnetic transition dipoles
				 dtensor3x3* ts, // quadrupole transition tensors
				 int number // number of interactions
				 )
// I put elements, dipoles, tensors  and wavevectors and orientation GetAt
{
	double	averaging=0;
	double 	bveraging[4]={0, 0, 0, 0};
	double 	mveraging[4]={0, 0, 0, 0};
	double dtemp;

	int matrix4[3][3] = {
		{ 4, -1, -1},
		{-1,  4, -1},
		{-1, -1,  4}
	};
	int matrix5[6][6] = {
		{ 3, -1, -1,  1,  1,  0},
		{-1,  3, -1, -1,  0,  1},
		{-1, -1,  3,  0, -1, -1},
		{ 1, -1,  0,  3, -1,  1},
		{ 1,  0, -1, -1,  3, -1},
		{ 0,  1, -1,  1, -1,  3}
	};
	double rvec[6];
	double lvec[6];

	// wavevectors and polarizations are directly specified.
	if(type == 0)// isotropic averaging
	{
		if(number == 4)// four vectors averaged
		{
			
			//////////////////////////////////////////////////////
			// DIPOLE APPROXIMATION
			
			// setting zeros;
			memset(rvec,0,6*sizeof(double));
			memset(lvec,0,6*sizeof(double));
					
			// Getting values:
			rvec[0] = (ds[3]*ds[2]) * (ds[1]*ds[0]);
			rvec[1] = (ds[3]*ds[1]) * (ds[2]*ds[0]);
			rvec[2] = (ds[3]*ds[0]) * (ds[2]*ds[1]);
			
			lvec[0] = (es[3]*es[2]) * (es[1]*es[0]);
			lvec[1] = (es[3]*es[1]) * (es[2]*es[0]);
			lvec[2] = (es[3]*es[0]) * (es[2]*es[1]);
			
					
			averaging=0;
			for(int ii=0;ii<3;ii++)
			for(int ik=0;ik<3;ik++)
				averaging += lvec[ii]*matrix4[ii][ik]*rvec[ik];
			
			averaging /= 30.0;


			if(ts) //  I go beyond the dipole approximation
			{
				for(int z=0; z<4; z++)
				{
					// this loop will select the wavevector term

				memset(rvec,0,6*sizeof(double));
				memset(lvec,0,6*sizeof(double));


				// Getting values:
				for(int l5=0;l5<3;l5++) // for 1-st expansion term
				for(int l4=0;l4<3;l4++) 
				for(int l3=0;l3<3;l3++)
				for(int l2=0;l2<3;l2++)
				for(int l1=0;l1<3;l1++)
//				for(int l0=0;l0<3;l0++)
				{
//					if(indvec == 3)
//						dtemp = ts[3].GetAt(l4*3+l3) *ds[2].GetAt(l2) *ds[1].GetAt(l1) *ds[0].GetAt(l0);
//					else if(indvec == 2)
//						dtemp = ds[3].GetAt(l3) *ts[2].GetAt(l4*3+l2) *ds[1].GetAt(l1) *ds[0].GetAt(l0);
//					else if(indvec == 1)
//						dtemp = ds[3].GetAt(l3) *ds[2].GetAt(l2) *ts[1].GetAt(l4*3+l1) *ds[0].GetAt(l0);
//					else if(indvec == 0)
//						dtemp = ds[3].GetAt(l3) *ds[2].GetAt(l2) *ds[1].GetAt(l1) *ts[0].GetAt(l4*3+l0);
//
//					if(l0==l4) rvec[0] += constants::lcs[l3][l2][l1]*dtemp;
//					if(l1==l4) rvec[1] += constants::lcs[l3][l2][l0]*dtemp;
//					if(l1==l0) rvec[2] += constants::lcs[l3][l2][l4]*dtemp;
//					if(l2==l4) rvec[3] += constants::lcs[l3][l1][l0]*dtemp;
//					if(l2==l0) rvec[4] += constants::lcs[l3][l1][l4]*dtemp;
//					if(l2==l1) rvec[5] += constants::lcs[l3][l0][l4]*dtemp;
//				}
//				lvec[0] = (es[3]*(es[2]^es[1])) * (es[0]*ks[indvec]);
//				lvec[1] = (es[3]*(es[2]^es[0])) * (es[1]*ks[indvec]);
//				lvec[2] = (es[3]*(es[2]^ks[indvec])) * (es[1]*es[0]);
//				lvec[3] = (es[3]*(es[1]^es[0])) * (es[2]*ks[indvec]);
//				lvec[4] = (es[3]*(es[1]^ks[indvec])) * (es[2]*es[0]);
//				lvec[5] = (es[3]*(es[0]^ks[indvec])) * (es[2]*es[1]);


					if(z==0) dtemp = constants::lcs[l5][l4][l3]*ds[3].GetAt(l4) *ds[2].GetAt(l3) *ds[1].GetAt(l2) *ts[0].GetAt(l5*3+l1);
					if(z==1) dtemp = constants::lcs[l5][l4][l3]*ds[3].GetAt(l4) *ds[2].GetAt(l3) *ts[1].GetAt(l5*3+l2) *ds[0].GetAt(l1);
					if(z==2) dtemp = constants::lcs[l5][l4][l3]*ds[3].GetAt(l4) *ts[2].GetAt(l5*3+l3) *ds[1].GetAt(l2) *ds[0].GetAt(l1);
					if(z==3) dtemp = constants::lcs[l5][l4][l3]*ts[3].GetAt(l5*3+l4) *ds[2].GetAt(l3) *ds[1].GetAt(l2) *ds[0].GetAt(l1);
					
					if(l2==l1) rvec[0] += constants::lcs[l5][l4][l3]*dtemp;
					if(l3==l1) rvec[1] += constants::lcs[l5][l4][l2]*dtemp;
					if(l3==l2) rvec[2] += constants::lcs[l5][l4][l1]*dtemp;
					if(l4==l1) rvec[3] += constants::lcs[l5][l3][l2]*dtemp;
					if(l4==l2) rvec[4] += constants::lcs[l5][l3][l1]*dtemp;
					if(l4==l3) rvec[5] += constants::lcs[l5][l2][l1]*dtemp;
				}
				
				lvec[0] = (ks[z]*(es[3]^es[2])) * (es[1]*es[0]);
				lvec[1] = (ks[z]*(es[3]^es[1])) * (es[2]*es[0]);
				lvec[2] = (ks[z]*(es[3]^es[0])) * (es[2]*es[1]);
				lvec[3] = (ks[z]*(es[2]^es[1])) * (es[3]*es[0]);
				lvec[4] = (ks[z]*(es[2]^es[0])) * (es[3]*es[1]);
				lvec[5] = (ks[z]*(es[1]^es[0])) * (es[3]*es[2]);

				bveraging[z]=0;
				for(int ii=0;ii<6;ii++)
				for(int ik=0;ik<6;ik++)
					bveraging[z] += lvec[ii]*matrix5[ii][ik]*rvec[ik];	
				bveraging[z] *= om[z]/30.0;
				//bveraging[z] = -bveraging[3];
				}// end of loop over four wavevectors	
					
			}// electric dipole stuff finished
				
				// now magnetic stuff:
			if(ms) // only when it is allocated
			{
					
					// vector K4:
				for(int indvec=0; indvec<4;indvec++)
				{
					// this loop will select the wavevector term
					
					// setting zeros;
				memset(rvec,0,6*sizeof(double));
				memset(lvec,0,6*sizeof(double));
					
				// Getting values:
				if(indvec == 0)
				{
					rvec[0] = (ds[3]*ds[2]) * (ds[1]*ms[0]);
					rvec[1] = (ds[3]*ds[1]) * (ds[2]*ms[0]);
					rvec[2] = (ds[3]*ms[0]) * (ds[2]*ds[1]);
					
					lvec[0]  = (es[3]*es[2]) * ((ks[0]^es[1])*es[0]);
					lvec[1]  = (es[3]*es[1]) * ((ks[0]^es[2])*es[0]);
					lvec[2]  = ((ks[0]^es[3])*es[0]) * (es[2]*es[1]);
				}
				else if(indvec == 1)
				{
					rvec[0] = (ds[3]*ds[2]) * (ms[1]*ds[0]);
					rvec[1] = (ds[3]*ms[1]) * (ds[2]*ds[0]);
					rvec[2] = (ds[3]*ds[0]) * (ds[2]*ms[1]);
					
					lvec[0]  = (es[3]*es[2]) * (es[1]*(ks[1]^es[0]));
					lvec[1]  = ((ks[1]^es[3])*es[1]) * (es[2]*es[0]);
					lvec[2]  = (es[3]*es[0]) * ((ks[1]^es[2])*es[1]);
				}
				else if(indvec == 2)
				{
					rvec[0] = (ds[3]*ms[2]) * (ds[1]*ds[0]);
					rvec[1] = (ds[3]*ds[1]) * (ms[2]*ds[0]);
					rvec[2] = (ds[3]*ds[0]) * (ms[2]*ds[1]);
					
					lvec[0]  = ((ks[2]^es[3])*es[2]) * (es[1]*es[0]);
					lvec[1]  = (es[3]*es[1]) * (es[2]*(ks[2]^es[0]));
					lvec[2]  = (es[3]*es[0]) * (es[2]*(ks[2]^es[1]));


				}
				else if(indvec == 3)
				{
					rvec[0] = (ms[3]*ds[2]) * (ds[1]*ds[0]);
					rvec[1] = (ms[3]*ds[1]) * (ds[2]*ds[0]);
					rvec[2] = (ms[3]*ds[0]) * (ds[2]*ds[1]);
					
					lvec[0]  = (es[3]*(ks[3]^es[2])) * (es[1]*es[0]);
					lvec[1]  = (es[3]*(ks[3]^es[1])) * (es[2]*es[0]);
					lvec[2]  = (es[3]*(ks[3]^es[0])) * (es[2]*es[1]);
				}

					
				mveraging[indvec]=0;
				for(int ii=0;ii<3;ii++)
				for(int ik=0;ik<3;ik++)
					mveraging[indvec] += lvec[ii]*matrix4[ii][ik]*rvec[ik];	
				mveraging[indvec] *= constants::sign(om[indvec])/(constants::c_cmfs*30.0);

				}// end of loop over four wavevectors	

		
			}// end of magnetic contributions

			return averaging + coni*( bveraging[3]+bveraging[2]+bveraging[1]+bveraging[0]+mveraging[3]+mveraging[2]+mveraging[1]+mveraging[0]);
//			return averaging + coni*( bveraging[3]+bveraging[2]+bveraging[1]+bveraging[0]-mveraging[3]-mveraging[2]-mveraging[1]-mveraging[0]);
		}// four vector isotropic averaging finished


		if(number == 2)// two vectors averaged still isotropic averaging
		{
			
			//////////////////////////////////////////////////////
			// DIPOLE APPROXIMATION
			
			averaging = (es[1]*es[0]) * (ds[1]*ds[0]);
			averaging /= 3.0;

			// FIRST ORDER IN WAVEVECTOR A: ELECTRIC CONTRIBUTIONS
			if(ts)// electric dipole stuff
			{					
				// setting zeros;
				bveraging[1] = 0;

				for(int l4=0;l4<3;l4++)// coordinate
				for(int l1=0;l1<3;l1++)
				for(int l0=0;l0<3;l0++)
				{
					dtemp = ts[1].GetAt(l4*3+l1) *ds[0].GetAt(l0);
					bveraging[1] += constants::lcs[l4][l1][l0]*dtemp;
				}
				bveraging[1] *= (ks[1]*(es[1]^es[0]));

				bveraging[1] *= om[1]/6.0;
//				bveraging[1] /= 6.0;
//				bveraging[1] = -bveraging[1];
					
				// FOR K1:
				// setting zeros;
				bveraging[0] = 0;

				for(int l4=0;l4<3;l4++)// coordinate
				for(int l1=0;l1<3;l1++)
				for(int l0=0;l0<3;l0++)
				{
					dtemp = ts[0].GetAt(l4*3+l0) *ds[1].GetAt(l1);
					bveraging[0] += constants::lcs[l4][l1][l0]*dtemp;
				}
				bveraging[0] *= (ks[0]*(es[1]^es[0]));

				bveraging[0] *= om[0]/6.0;
//				bveraging[0] /= 6.0;
			}// electric dipole stuff finished
				
			// now magnetic stuff:
			if(ms) // only when it is allocated
			{
				// vector K2:
			
				mveraging[1] = (es[1]*(ks[1]^es[0])) * (ms[1]*ds[0]);
//				mveraging[1] /= 3.0;
//				mveraging[1] /= (constants::c_cmfs * ks[1].Amplitude() );
//				mveraging[1] = -mveraging[1];
				mveraging[1] *= constants::sign(om[1])/(constants::c_cmfs*3.0);
					
				// 
				// vector K1:
				
				mveraging[0] = ((ks[0]^es[1])*es[0]) * (ms[0]*ds[1]);
				mveraging[0] *= constants::sign(om[0])/(constants::c_cmfs*3.0);
//				mveraging[0] /= 3.0;
//				mveraging[0] /= (constants::c_cmfs * ks[0].Amplitude() );
//				mveraging[0] = -mveraging[0];
					
			}// end of magnetic contributions


//			return averaging + coni*( bveraging[1]+bveraging[0]-mveraging[1]-mveraging[0]);
			return averaging + coni*( bveraging[1]+bveraging[0]+mveraging[1]+mveraging[0]);
		}// two vector isotropic averaging finished


		if(number == 3)// three vectors averaged
		{
			//////////////////////////////////////////////////////
			// DIPOLE APPROXIMATION
			
			// setting zeros;
			averaging = (es[2]*(es[1]^es[0])) * (ds[2]*(ds[1]^ds[0]));
			averaging /= 6.0;
			return averaging;

			// magnetic dipoles, etc will not be used
			
		}// three vector isotropic averaging finished
		
						
	}// isotropic parts finished
	else if(type == 1)// explicit orientation
	{
		if(number == 4)// four vectors averaged
		{
			return (ds[3]*es[3])*(ds[2]*es[2])*(ds[1]*es[1])*(ds[0]*es[0]);
		}
		else if(number == 3)
		{
			return (ds[2]*es[2])*(ds[1]*es[1])*(ds[0]*es[0]);
		}
		else if(number == 2)
		{
			return (ds[1]*es[1])*(ds[0]*es[0]);
		}
		else return 0.0;

	}
	else 
		std::cout<<"Error: this orientational averaging is not coded.\n";

	return 0.0;

}



