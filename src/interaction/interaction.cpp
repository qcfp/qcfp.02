#include"interaction.hpp"

interaction::interaction()
{
	edipole = new dvector3d[1];
	efield = new dvector3d[1];
	kvector = new dvector3d[1];
	omega = new double[1];

	approach = 0;
	if(approach)
	{
		mdipole = new dvector3d[1];
		quadrupole = new dtensor3x3[1];
	}
	size = 1;
	aamplitude =0.0;
	aprocess.SetType(0);

	trl = new int[1];
	trr = new int[1];
	reassignment = 1;

	sysd = 0;
	syst = 0;
	sysm = 0;
}

interaction::interaction(int length, int level,int avtype)
{
	if(length <=0)
	{
		std::cout<<"Error: incorrect parameter in interaction::size\n";
		return;
	}
	size = length;
	edipole = new dvector3d[size];
	efield = new dvector3d[size];
	kvector = new dvector3d[size];
	omega = new double[size];

	

	approach = level;
	if(approach)
	{
		mdipole = new dvector3d[size];
		quadrupole = new dtensor3x3[size];
	}

	aamplitude = 0.0;
	aprocess.SetType(avtype);

	trl = new int[size];
	trr = new int[size];
	reassignment = 1;

	sysd = 0;
	syst = 0;
	sysm = 0;
}

interaction::~interaction()
{
	delete[] edipole;// = new dvector3d[1];
	delete[] efield;// = new dvector3d[1];
	delete[] kvector;// = new dvector3d[1];
	delete[] omega;// = new double[1];

	if(approach)
	{
		delete[] mdipole;// = new dvector3d[1];
		delete[] quadrupole;// = new dtensor3x3[1];
	}

	delete[] trl;
	delete[] trr;

	sysd = 0;
	syst = 0;
	sysm = 0;

}
complexv interaction::GetAveragedAmplitude()
{
	// calculation should proceed here
	// 
	dtensor3x3* loct=0;
	dvector3d* locm=0;
	if(reassignment)
	{
		if(syst) loct = quadrupole;
		if(sysm) locm = mdipole;
		for(int iii=0; iii<size; iii++){
		edipole[iii]=sysd[trl[iii]][trr[iii]];
		if(syst) loct[iii]=syst[trl[iii]][trr[iii]];
		if(sysm) locm[iii]=sysm[trl[iii]][trr[iii]];
		}
		reassignment = 0;

		aamplitude = aprocess.rot_av_complete_2_4(
			efield, // optical polarizations
			kvector, // optical wavevectors
			omega,  // fourier frequencies
			edipole, // electric transition dipoles
			locm, // magnetic transition dipoles
			loct, // quadrupole transition tensors
			size // number of interactions
			 );
	}
	return aamplitude;

};

void interaction::PopulateSys(dvector3d** isysd,dtensor3x3** isyst,dvector3d** isysm)
{
	sysd=isysd;
	syst=isyst;
	sysm=isysm;
}

void interaction::PopulateFields(dvector3d* ie,dvector3d* ik,double* io)
{
	for(int ind=0; ind<size; ind++)
	{
		efield[ind]=ie[ind];
		kvector[ind]=ik[ind];
		omega[ind]=io[ind];
	}
}

void interaction::AssignLevels(int number,int itl,int itr)
{
	trl[number]=itl;
	trr[number]=itr;

	reassignment = 1;
}

