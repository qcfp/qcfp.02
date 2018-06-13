#include"constants.hpp"


double constants::pi = 3.141592653589793;
double constants::pi2 = 2.0*constants::pi;
double constants::piover2 = constants::pi/2.0;
double constants::root2 = 1.4142135623730951;
double constants::root3 = 1.7320508075688772;
double constants::cos45 = 0.707106781186548;

double constants::smallepsilon = 1e-13;

int constants::lcs2[2][2] = {
		{ 0, 1},
		{-1,  0} };
int constants::lcs[3][3][3] = {
		{{0, 0, 0},
		 {0, 0, 1},
		 {0, -1, 0}},
		{{0, 0, -1},
		 {0, 0, 0},
		 {1, 0, 0}},
		{{0, 1, 0},
		 {-1, 0, 0},
		 {0, 0, 0}}  };



// physical constants
double constants::c_cmfs = 2.99792458e-5; // speed of light [cm/fs]
double constants::c_cmfs_r = constants::c_cmfs*constants::pi2; // radial speed of light 2*pi*c [cm/fs] 
double constants::c_cmfs_i = 1.0/constants::c_cmfs_r; // inverse radial speed of light 1/(2*pi*c) [fs/cm] 
double constants::bkTc = 0.695034766363; //boltzmann constant [cm-1/K]
double constants::bkTe = 8.61733247878e-5; //boltzmann constant [eV/K]

