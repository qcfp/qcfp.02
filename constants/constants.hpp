#pragma once

class constants
{
public:
    
// mathematical constants
     static double pi;
     static double piover2;
     static double pi2;
     static double root2;
     static double root3;
     static double cos45;

	// levi civita symbol
	static int lcs2[2][2];
	static int lcs[3][3][3];

    static double smallepsilon;
    

    // physical constants
     static double c_cmfs; // speed of light [cm/fs]
     static double c_cmfs_r; // radial speed of light 2*pi*c [cm/fs] 
     static double c_cmfs_i; // inverse radial speed of light 1/(2*pi*c) [fs/cm] 
     static double bkTc; //boltzmann constant [cm-1/K]
     static double bkTe; //boltzmann constant [eV/K]

    
/*    constants(){
        pi = 3.141592653589793;
        pi2 = 2.0*pi;
        root2 = 1.4142135623730951;
        root3 = 1.7320508075688772;
        
        // special coefficient cos(45)
        cos45 = 0.707106781186548;
        
        c_cmfs = 2.99792458e-5;
        c_cmfs_r = c_cmfs*pi2;
        c_cmfs_i = 1.0/c_cmfs_r;
        bkTc = 0.695034766363;
        bkTe = 8.61733247878e-5;
    }
 */
    
	public: inline static double sign(double x)
	{
		return (x>0.0) - (x<0.0);
	}
};


