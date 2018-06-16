#include<ios>
#include"../include/conversions.h"
#include"../include/universe.hpp"



r_tensor3x3 UnitsConversion(
	spectron_units units, 
	r_tensor3x3  vector,
	spectron_property prop, 
	double parameter
	)
{
	if(prop == prop_quad)
	{
		if( units == units_local ) //1
			// on the input must be [A2*e] ; no conversion 
			return vector;

		if( units == units_de_nm_cm ) //2
			// on the input must be [Debye nm]. Converting to [A2*e] 
			return vector*2.0819435;

		if( units == units_de_an_cm ) //3
			// on the input must be [Debye A]. Converting to [A2*e] 
			return vector*0.20819435;

		if( units == units_bm_nm_ev ) //4
			// on the input must be [BohrMagneton(c) nm]. Converting to [A2*e] 
			return vector*0.019308;

		if( units == units_bm_an_ev ) //5
			// on the input must be [BohrMagneton(c) A]. Converting to [A2*e] 
			return vector*0.0019308;

		if( units == units_db_an_cm ) //6
			// on the input must be [Debye A]. Converting to [A2*e] 
			return vector*0.20819435;

		if( units == units_db_an_ev ) //7
			// on the input must be [Debye A]. Converting to [A2*e] 
			return vector*0.20819435;

		if( units == units_db_nm_ev ) //8
			// on the input must be [Debye nm]. Converting to [A2*e] 
			return vector*2.0819435;

		if( units == units_ia_an_ev ) //9
			// on the input must be [A2*e] ; no conversion 
			return vector;

		if( units == units_ia_an_cm ) //10
			// on the input must be [A2*e] ; no conversion 
			return vector;

		if( units == units_de_nm_th ) //11
			// on the input must be [Debye nm]. Converting to [A2*e] 
			return vector*2.0819435;

		if( units == units_bm_nm_th ) //12
			// on the input must be [BohrMagneton(c) nm]. Converting to [A2*e] 
			return vector*0.019308;
	}
}


r_vector3D UnitsConversion(
	spectron_units units, 
	r_vector3D vector,
	spectron_property prop, 
	double parameter
	)
{
	if(prop == prop_edip)
	{
		if( units == units_local ) //1
			// R*e representation: on the input must be [A*e] ; no conversion 
			return vector;

		if( units == units_de_nm_cm ) //2
			// R*e representation: on the input must be [Debye]. Converting to [A*e] 
			return vector*0.20819435;

		if( units == units_de_an_cm ) //3
			// R*e representation: on the input must be [Debye]. Converting to [A*e] 
			return vector*0.20819435;

		if( units == units_bm_nm_ev ) //4
			// R*e representation: on the input must be [BohrMagneton(c)]. Converting to [A*e] 
			return vector*0.0019308;

		if( units == units_bm_an_ev ) //5
			// R*e representation: on the input must be [BohrMagneton(c)]. Converting to [A*e] 
			return vector*0.0019308;

		if( units == units_db_an_cm ) //6
			// R*e representation: on the input must be [Debye]. Converting to [A*e] 
			return vector*0.20819435;

		if( units == units_db_an_ev ) //7
			// R*e representation: on the input must be [Debye]. Converting to [A*e] 
			return vector*0.20819435;

		if( units == units_db_nm_ev ) //8
			// R*e representation: on the input must be [Debye]. Converting to [A*e] 
			return vector*0.20819435;

		if( units == units_ia_an_ev ) //9
			// nabla representation: on the input must be [A^{-1}]. Converting to [A*e] 
			return vector*( 6.1459e4 / parameter);

		if( units == units_ia_an_cm ) //10
			// nabla representation: on the input must be [A^{-1}]. Converting to [A*e] 
			return vector*( 6.1459e4 / parameter);

		if( units == units_de_nm_th ) //11
			// R*e representation: on the input must be [Debye]. Converting to [A*e] 
			return vector*0.20819435;

		if( units == units_bm_nm_th ) //12
			// R*e representation: on the input must be [BohrMagneton(c)]. Converting to [A*e] 
			return vector*0.0019308;

	}
	else if(prop == prop_coor)
	{

		if( units == units_local ) //1
			// on the input must be [A] ; no conversion 
			return vector;

		if( units == units_de_nm_cm ) //2
			// on the input must be [nm]. Converting to [A] 
			return vector*10.0;

		if( units == units_de_an_cm ) //3
			// on the input must be [A] ; no conversion 
			return vector;

		if( units == units_bm_nm_ev ) //4
			// on the input must be [nm]. Converting to [A] 
			return vector*10.0;

		if( units == units_bm_an_ev ) //5
			// on the input must be [A] ; no conversion 
			return vector;

		if( units == units_db_an_cm ) //6
			// on the input must be [A] ; no conversion 
			return vector;

		if( units == units_db_an_ev ) //7
			// on the input must be [A] ; no conversion 
			return vector;

		if( units == units_db_nm_ev ) //8
			// on the input must be [nm]. Converting to [A] 
			return vector*10.0;

		if( units == units_ia_an_ev ) //9
			// on the input must be [A] ; no conversion 
			return vector;

		if( units == units_ia_an_cm ) //10
			// on the input must be [A] ; no conversion 
			return vector;

		if( units == units_de_nm_th ) //11
			// on the input must be [nm]. Converting to [A] 
			return vector*10.0;

		if( units == units_bm_nm_th ) //12
			// on the input must be [nm]. Converting to [A] 
			return vector*10.0;

	}
	else if(prop == prop_mdip)
	{

		if( units == units_local ) //1
			// on the input must be [A^2*e/fs] 
			return vector;

		if( units == units_de_nm_cm ) //2
			// on the input must be [Debye] (due to speed of light). Converting to [A^2*e/fs] (no speed of light)
			return vector*(624.16);

		if( units == units_de_an_cm ) //3
			// on the input must be [Debye] (due to speed of light). Converting to [A^2*e/fs] (no speed of light)
			return vector*(624.16);

		if( units == units_bm_nm_ev ) //4
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

		if( units == units_bm_an_ev ) //5
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

		if( units == units_db_an_cm ) //6
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

		if( units == units_db_an_ev ) //7
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

		if( units == units_db_nm_ev ) //8
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

		if( units == units_ia_an_ev ) //9
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

		if( units == units_ia_an_cm ) //10
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

		if( units == units_de_nm_th ) //11
			// on the input must be [Debye] (due to speed of light). Converting to [A^2*e/fs] (no speed of light)
			return vector*(624.16);

		if( units == units_bm_nm_th ) //12
			// on the input must be [BohrMagneton]. Converting to [A^2*e/fs] (no speed of light)
			return vector*(5.7883);

	}

}



double UnitsConversion(
	spectron_units units, 
	double value,
	spectron_property prop, 
	double parameter
	)
{
	if(prop == prop_ener)// energies
	{

		if( units == units_local ) //1
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_de_nm_cm ) //2
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_de_an_cm ) //3
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_bm_nm_ev ) //4
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value*8065.545;

		if( units == units_bm_an_ev ) //5
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value*8065.545;

		if( units == units_db_an_cm ) //6
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_db_an_ev ) //7
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value*8065.545;

		if( units == units_db_nm_ev ) //8
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value*8065.545;

		if( units == units_ia_an_ev ) //9
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value*8065.545;

		if( units == units_ia_an_cm ) //10
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_de_nm_th ) //11
			// on the input must be [THz]. converting to [cm^{-1}]  
			return value*33.357;

		if( units == units_bm_nm_th ) //12
			// on the input must be [THz]. converting to [cm^{-1}]  
			return value*33.357;
	}
	if(prop == prop_dist)// distances
	{

		if( units == units_local ) //1
			// on the input must be A ; no conversion 
			return value;

		if( units == units_de_nm_cm ) //2
			// on the input must be nm ; converting to A 
			return value*10.0;

		if( units == units_de_an_cm ) //3
			// on the input must be A ; no conversion 
			return value;

		if( units == units_bm_nm_ev ) //4
			// on the input must be nm ; converting to A 
			return value*10.0;

		if( units == units_bm_an_ev ) //5
			// on the input must be A ; no conversion 
			return value;

		if( units == units_db_an_cm ) //6
			// on the input must be A ; no conversion 
			return value;

		if( units == units_db_an_ev ) //7
			// on the input must be A ; no conversion 
			return value;

		if( units == units_db_nm_ev ) //8
			// on the input must be nm ; converting to A 
			return value*10.0;

		if( units == units_ia_an_ev ) //9
			// on the input must be A ; no conversion 
			return value;

		if( units == units_ia_an_cm ) //10
			// on the input must be A ; no conversion 
			return value;

		if( units == units_de_nm_th ) //11
			// on the input must be nm ; converting to A 
			return value*10.0;

		if( units == units_bm_nm_th ) //12
			// on the input must be nm ; converting to A 
			return value*10.0;

	}
}

double UnitsConversion_Inv(
	spectron_units units, 
	double value,
	spectron_property prop, 
	double parameter
	)
{
	// all are inverse
	if(prop == prop_ener)// energies
	{

		if( units == units_local ) //1
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_de_nm_cm ) //2
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_de_an_cm ) //3
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_bm_nm_ev ) //4
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value/8065.545;

		if( units == units_bm_an_ev ) //5
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value/8065.545;

		if( units == units_db_an_cm ) //6
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_db_an_ev ) //7
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value/8065.545;

		if( units == units_db_nm_ev ) //8
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value/8065.545;

		if( units == units_ia_an_ev ) //9
			// on the input must be [eV]. converting to [cm^{-1}]  
			return value/8065.545;

		if( units == units_ia_an_cm ) //10
			// on the input must be [cm^{-1}] ; no conversion 
			return value;

		if( units == units_de_nm_th ) //11
			// on the input must be [THz]. converting to [cm^{-1}]  
			return value/33.357;

		if( units == units_bm_nm_th ) //12
			// on the input must be [THz]. converting to [cm^{-1}]  
			return value/33.357;
	}
	if(prop == prop_dist)// distances
	{
		if( units == units_local ) //1
			// on the input must be A ; no conversion 
			return value;

		if( units == units_de_nm_cm ) //2
			// on the input must be nm ; converting to A 
			return value/10.0;

		if( units == units_de_an_cm ) //3
			// on the input must be A ; no conversion 
			return value;

		if( units == units_bm_nm_ev ) //4
			// on the input must be nm ; converting to A 
			return value/10.0;

		if( units == units_bm_an_ev ) //5
			// on the input must be A ; no conversion 
			return value;

		if( units == units_db_an_cm ) //6
			// on the input must be A ; no conversion 
			return value;

		if( units == units_db_an_ev ) //7
			// on the input must be A ; no conversion 
			return value;

		if( units == units_db_nm_ev ) //8
			// on the input must be nm ; converting to A 
			return value/10.0;

		if( units == units_ia_an_ev ) //9
			// on the input must be A ; no conversion 
			return value;

		if( units == units_ia_an_cm ) //10
			// on the input must be A ; no conversion 
			return value;

		if( units == units_de_nm_th ) //11
			// on the input must be nm ; converting to A 
			return value/10.0;

		if( units == units_bm_nm_th ) //12
			// on the input must be nm ; converting to A 
			return value/10.0;

	}
}


////////////////////////////////
r_vector3D UnitsConversion(
	spectron_units units, 
	r_vector3D vector,
	spectron_property prop 
	)
{
	return 
	UnitsConversion(units,vector,prop,0.0);
}

r_tensor3x3 UnitsConversion(
	spectron_units units, 
	r_tensor3x3 vector,
	spectron_property prop 
	)
{
	return 
	UnitsConversion(units,vector,prop,0.0);
}

double UnitsConversion(
	spectron_units units, 
	double value,
	spectron_property prop 
	)
{
	return 
	UnitsConversion(units,value,prop,0.0);
}

double UnitsConversion_Inv(
	spectron_units units, 
	double value,
	spectron_property prop
	)
{
	return 
	UnitsConversion_Inv(units,value,prop,0.0);
}




