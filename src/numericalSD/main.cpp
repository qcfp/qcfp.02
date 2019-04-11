// doing some preliminary testings

#include"../storage/storage.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../interpolationF/interpolationF.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../complexv/complexv.hpp"
#include"../toolsFFT/toolsFFT.hpp"
#include"../constants/constants.hpp"

#include"numericalSD.hpp"




#include<string>
using std::string;


int main(int argc, const char * argv[])
{

    if(argc != 3 )
    {    cout<<"Error: specify the spectral density file name and temperature\n";
        return 0;
    }
 
    toolsIO tio;
    
    string filename(argv[1]);
    double temperature = tio.fromString<double>(string(argv[2]));
    
	// creating numeric SD object
	numericalSD* obj = new numericalSD(temperature);
	//numericalSD* obj = new numericalSD(77);
	// that call does not do any precalculations
	// it only reads the spectral density from a file
    //obj->updateUnits("cm-fs-K"); // su units negerai veikia
    obj->ReadSD(filename);


    
    asymptoticLF_complexv gfun = obj->GetGf();
    asymptoticLF_complexv mfun = obj->GetMf();
    asymptoticLF_complexv cfun = obj->GetCf();

    
    obj->printACCESSORY(filename+".acce.txt");
    
    
    gfun.SaveF(filename+".gfun.wrk");
    mfun.SaveF(filename+".mfun.wrk");
    cfun.SaveF(filename+".cfun.wrk");
    
    //	// returning lineshape function
//	comple gfun = obj->ReturnLSFunct(0.0);

	// writing out M function
	obj->printMBAR(filename+".mfun.txt");

	// writing out G function
	obj->printGFUN(filename+".gfun.txt");

	// writing out Access functions
	obj->printACCESSORY(filename+".acce.txt");

	// writing out Correlation function
	obj->printCORELF(filename+".cfun.txt");

        
//        // making binary outputs for high precision I/O
//        asymptoticLF_complexv cfun;
//        asymptoticLF<double> rfun;
//        ofstream ofs;
//        string fnb;
//        
//        // spectral density
//        rfun = obj->GetSD();
//        fnb = filename+"-Spde.bin";
//        ofs.open(fnb.c_str());
//        rfun.SaveF(&ofs);
//        ofs.close();
//        
//        // correlation function
//        cfun = obj->GetCf();
//        fnb = filename+"-Corr.bin";
//        ofs.open(fnb.c_str());
//        cfun.SaveF(&ofs);
//        ofs.close();
//        
//        // lineshape function
//        cfun = obj->GetGf();
//        fnb = filename+"-Gfun.bin";
//        ofs.open(fnb.c_str());
//        cfun.SaveF(&ofs);
//        ofs.close();
    
        
        
	// cleaning
	delete obj;

}

