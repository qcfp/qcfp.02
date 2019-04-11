
#ifndef _CONVERSIONS_H_
#define _CONVERSIONS_H_


#include "universe.hpp"
#include "constants.h"
#include "../../libshared/include/simple-vectors.h"
#include "../../libshared/include/simple-tensors.h"



//#define HAR_2_CM1 219474.63
//#define CM1_2_HAR 1 / HAR_2_CM1

//#define HAR_2_EV 27.2113961
//#define EV_2_HAR 1 / HAR_2_EV

//#define HAR_2_KCPM 627.510
//#define KCPM_2_HAR 1 / HAR_2_KCPM

//#define HAR_2_K 315773
//#define K_2_HAR 1 / HAR_2_K




r_vector3D UnitsConversion(
	spectron_units units, 
	r_vector3D vector,
	spectron_property prop, 
	double parameter
	);

r_tensor3x3 UnitsConversion(
	spectron_units units, 
	r_tensor3x3 vector,
	spectron_property prop, 
	double parameter
	);

double UnitsConversion(
	spectron_units units, 
	double value,
	spectron_property prop, 
	double parameter
	);

double UnitsConversion_Inv(
	spectron_units units, 
	double value,
	spectron_property prop, 
	double parameter
	);
////////////////////////////////
r_vector3D UnitsConversion(
	spectron_units units, 
	r_vector3D vector,
	spectron_property prop 
	);

r_tensor3x3 UnitsConversion(
	spectron_units units, 
	r_tensor3x3 vector,
	spectron_property prop
	);

double UnitsConversion(
	spectron_units units, 
	double value,
	spectron_property prop 
	);

double UnitsConversion_Inv(
	spectron_units units, 
	double value,
	spectron_property prop
	);




#endif // _CONVERSIONS_H_
