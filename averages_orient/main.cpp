// example

#include<iostream>

#include"averages_orient.hpp"
#include"../toolsRandom/toolsRandom.hpp"

// this example program shows how one gets the isotropically averaged
// dipole products

// local frame unit dipoles are generated randomly
// lab frame dipoles cover the possible independent configurations

using std::cin;
using std::cout;

int main(int argc, const char * argv[])
{
    toolsRandom trd;
    averages_orient aor(0);
    dvector3d vec[4];
    dvector3d eve[4];

	vec[0] = trd.RandV3D();
	vec[1] = trd.RandV3D();
	vec[2] = trd.RandV3D();
	vec[3] = trd.RandV3D();
    


	cout<<"creating two-vector isotropic orientational averaging\n";

	cout<<"1a. Making two unit-length random vectors:\n";
	cout<<"vec1: "<<vec[0].x()<<" "<<vec[0].y()<<" "<<vec[0].z()<<"\n";
	cout<<"vec2: "<<vec[1].x()<<" "<<vec[1].y()<<" "<<vec[1].z()<<"\n";

	cout<<"1b. their orientational averaging with X lab vector:\n";
	eve[0] = dvector3d(0.,0.,1.);
	eve[1] = dvector3d(0.,0.,1.);
	cout<<aor.rot_av_dip_2(eve,vec)<<"\n";


	cout<<"2a. Making four unit-length random vectors:\n";
	cout<<"vec1: "<<vec[0].x()<<" "<<vec[0].y()<<" "<<vec[0].z()<<"\n";
	cout<<"vec2: "<<vec[1].x()<<" "<<vec[1].y()<<" "<<vec[1].z()<<"\n";
	cout<<"vec3: "<<vec[2].x()<<" "<<vec[2].y()<<" "<<vec[2].z()<<"\n";
	cout<<"vec4: "<<vec[3].x()<<" "<<vec[3].y()<<" "<<vec[3].z()<<"\n";

	cout<<"2b. their orientational averaging with XXYY lab vectors:\n";
	eve[0] = dvector3d(1.,0.,0.);
	eve[1] = dvector3d(1.,0.,0.);
	eve[2] = dvector3d(0.,1.,0.);
	eve[3] = dvector3d(0.,1.,0.);
	cout<<aor.rot_av_dip_4(eve,vec)<<"\n";

	cout<<"2c. their orientational averaging with XYYX lab vectors:\n";
	eve[0] = dvector3d(1.,0.,0.);
	eve[1] = dvector3d(0.,1.,0.);
	eve[2] = dvector3d(0.,1.,0.);
	eve[3] = dvector3d(1.,0.,0.);
	cout<<aor.rot_av_dip_4(eve,vec)<<"\n";

	cout<<"2d. their orientational averaging with XYXY lab vectors:\n";
	eve[0] = dvector3d(1.,0.,0.);
	eve[1] = dvector3d(0.,1.,0.);
	eve[2] = dvector3d(1.,0.,0.);
	eve[3] = dvector3d(0.,1.,0.);
	cout<<aor.rot_av_dip_4(eve,vec)<<"\n";

	cout<<"2e. their orientational averaging with ZZZZ lab vectors:\n";
	eve[0] = dvector3d(0.,0.,1.);
	eve[1] = dvector3d(0.,0.,1.);
	eve[2] = dvector3d(0.,0.,1.);
	eve[3] = dvector3d(0.,0.,1.);
	cout<<aor.rot_av_dip_4(eve,vec)<<"\n";

	cout<<"That is all.\n";

    return 0;
}
