#pragma once 
#include"../dvector3d/dvector3d.hpp"
#include"../dtensor3x3/dtensor3x3.hpp"
#include"../complexv/complexv.hpp"
#include"../averages_orient/averages_orient.hpp"

class interaction
{
private:
	dvector3d* edipole;
	dvector3d* mdipole;
	dtensor3x3* quadrupole;

	dvector3d* efield;
	dvector3d* kvector;
	double* omega;

private:
	int size;
	int approach;
	complexv aamplitude;
	averages_orient aprocess;
	int* trl;
	int* trr;
	int reassignment;

	dvector3d** sysd;
	dvector3d** sysm;
	dtensor3x3** syst;

	// functions
public:
	interaction();
	~interaction();
	interaction(int length, int level,int avtype);

	void PopulateSys(dvector3d** isysd, dtensor3x3** isyst, dvector3d** isysm);
	void PopulateFields(dvector3d* ie,dvector3d* ik,double* io);
	void AssignLevels(int number,int itl,int itr);

	complexv GetAveragedAmplitude();
};



