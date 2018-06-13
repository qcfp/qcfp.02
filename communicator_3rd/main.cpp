//


#include"../storage/storage.hpp"
#include"../toolsIO/toolsIO.hpp"
#include "communicator_3rd.hpp"

#include<string>
using namespace std;
//  main.cpp
// empty
int main(int argc, const char * argv[])
{
    // here we will have means to convert excitonic model to the input of levels
    // bath will create wrk files so that later these can be used

    // all necessary information is read from a single file
    if (argc!=3){
            cout<<"Error: specify the input file and the output file\n";
            return 0;
    }
    string ifilename(argv[1]);
    string ofilename(argv[2]);

    
    
    ifstream istr(ifilename.c_str());
    if(!istr.is_open()){
        cout<<"Error: communicator_3rd file has not been opened\n";
        return 0;
    }

    
    toolsIO tio;
    communicator_3rd calculator;
    calculator.eigensys=false;
	// reading the whole input file
	tio.ReadWholeFile(ifilename,calculator.inputfile,1);
    tio.StreamSkipTrailers(istr);
    string tstr; getline(istr,tstr);
    calculator.ReadSystem(istr);
    calculator.ReadBath(istr);

    calculator.MakeESSystem();


    // writing the resulting system to the output file
    std::ofstream fo(ofilename.c_str());
    
    fo.precision(14);
    
    fo<<"#########################################\n";
    fo<<calculator.complexity<<"\n";
    fo<<calculator.bath_complexity<<"\n";
    
    
    if(calculator.band == 2){
    fo<<"EigensysNumberOf012Levels:\n";
    fo<<"# number of zero-band levels:\n";
    fo<<1<<"\n";
    fo<<"# number of one-band levels:\n";
    fo<<calculator.numE<<"\n";
    fo<<"# number of two-band levels:\n";
    fo<<calculator.numF<<"\n";
    } else {
    fo<<"EigensysNumberOf01Levels:\n";
    fo<<"# number of zero-band levels:\n";
    fo<<1<<"\n";
    fo<<"# number of one-band levels:\n";
    fo<<calculator.numE<<"\n";
    }
    // print eigenvalues
    fo<<"EigensysLevelEnergies:\n";
    fo<<"# zero-band level energies:\n";
    fo<<0.0<<"\n";

    // print eigenvalues
    fo<<"# one-band level energies:\n";
    for(int ind=0; ind<calculator.numE; ind++)
        fo<<calculator.evals.data1D[ind+1]<<"\n";

    if(calculator.numF>0){
    // print eigenvalues 2 
    fo<<"# two-band level energies:\n";
    for(int ind=0; ind<calculator.numF; ind++)
        fo<<calculator.evals.data1D[ind+calculator.numE+1]<<"\n";
    }

    
    
    
    fo<<"# Transition dipoles:\n";
    fo<<"EigensysTransitionDipoles:\n";
    fo<<"# number of entries:\n";
    fo<<calculator.numE+calculator.numE*calculator.numF<<"\n";


    // print eigendipoles
    fo<<"# dipoles G-E:\n";
    for(int ind=0; ind<calculator.numE; ind++)
    {
        fo<<"0\t"<<ind+1<<"\t";
        fo<<calculator.edips.data2D[0][ind+1].x()<<"\t";
        fo<<calculator.edips.data2D[0][ind+1].y()<<"\t";
        fo<<calculator.edips.data2D[0][ind+1].z()<<"\n";
    }

    // print eigendipoles
    if(calculator.numF>0){
    fo<<"# dipoles E-F:\n";
    for(int ind=0; ind<calculator.numE; ind++)
        for(int inf=0; inf<calculator.numF; inf++)
        {
            fo<<ind+1<<"\t"<<inf+1+calculator.numE<<"\t";
            //        cout<<ind<<"\t"<<inf<<"\t";
            fo<<calculator.edips.data2D[ind+1][inf+calculator.numE+1].x()<<"\t";
            fo<<calculator.edips.data2D[ind+1][inf+calculator.numE+1].y()<<"\t";
            fo<<calculator.edips.data2D[ind+1][inf+calculator.numE+1].z()<<"\n";
        }
    }

    // next magnetic dipoles
    if(calculator.ed_m.IsAlloc()){

    fo<<"# Magnetic transition dipoles:\n";
    fo<<"EigensysMagneticDipoles:\n";
    fo<<"# number of entries:\n";
    fo<<calculator.numE+calculator.numE*calculator.numF<<"\n";


    // print eigendipoles
    fo<<"# dipoles G-E:\n";
    for(int ind=0; ind<calculator.numE; ind++)
    {
        fo<<"0\t"<<ind+1<<"\t";
        fo<<calculator.ed_m.data2D[0][ind+1].x()<<"\t";
        fo<<calculator.ed_m.data2D[0][ind+1].y()<<"\t";
        fo<<calculator.ed_m.data2D[0][ind+1].z()<<"\n";
    }

    // print eigendipoles
    if(calculator.numF>0){
    fo<<"# dipoles E-F:\n";
    for(int ind=0; ind<calculator.numE; ind++)
        for(int inf=0; inf<calculator.numF; inf++)
        {
            fo<<ind+1<<"\t"<<inf+1+calculator.numE<<"\t";
            //        cout<<ind<<"\t"<<inf<<"\t";
            fo<<calculator.ed_m.data2D[ind+1][inf+calculator.numE+1].x()<<"\t";
            fo<<calculator.ed_m.data2D[ind+1][inf+calculator.numE+1].y()<<"\t";
            fo<<calculator.ed_m.data2D[ind+1][inf+calculator.numE+1].z()<<"\n";
        }
    }
    }
    
    // next quadrupoles dipoles
    if(calculator.eten.IsAlloc()){

    fo<<"# Electric quadrupoles of transitions:\n";
    fo<<"EigensysQuadrupoles:\n";
    fo<<"# number of entries:\n";
    fo<<calculator.numE+calculator.numE*calculator.numF<<"\n";


    // print eigen quadrupoles
    fo<<"# quadrupoles (xx,xy,xz,yy,yz,zz) G-E:\n";
    for(int ind=0; ind<calculator.numE; ind++)
    {
        fo<<"0\t"<<ind+1<<"\t";
        fo<<calculator.eten.data2D[0][ind+1].xx()<<"\t";
        fo<<calculator.eten.data2D[0][ind+1].xy()<<"\t";
        fo<<calculator.eten.data2D[0][ind+1].xz()<<"\t";
        fo<<calculator.eten.data2D[0][ind+1].yy()<<"\t";
        fo<<calculator.eten.data2D[0][ind+1].yz()<<"\t";
        fo<<calculator.eten.data2D[0][ind+1].zz()<<"\n";
    }

    // print eigendipoles
    if(calculator.numF>0){
    fo<<"# quadrupoles E-F:\n";
    for(int ind=0; ind<calculator.numE; ind++)
        for(int inf=0; inf<calculator.numF; inf++)
        {
            fo<<ind+1<<"\t"<<inf+1+calculator.numE<<"\t";
            //        cout<<ind<<"\t"<<inf<<"\t";
            fo<<calculator.eten.data2D[ind+1][inf+calculator.numE+1].xx()<<"\t";
            fo<<calculator.eten.data2D[ind+1][inf+calculator.numE+1].xy()<<"\t";
            fo<<calculator.eten.data2D[ind+1][inf+calculator.numE+1].xz()<<"\t";
            fo<<calculator.eten.data2D[ind+1][inf+calculator.numE+1].yy()<<"\t";
            fo<<calculator.eten.data2D[ind+1][inf+calculator.numE+1].yz()<<"\t";
            fo<<calculator.eten.data2D[ind+1][inf+calculator.numE+1].zz()<<"\n";
        }
    }
    }
    
    
    
    
    fo<<"BathTemperature:\n";
    fo<<"# Defining the bath and SB-couplings\n";
    fo<<"# Temperature\n";
    fo<<calculator.tempr<<"\n";
    
    //if(calculator.complexity.compare("key-setup-complete")==0)
    //{
    //    fo<<"# ZPL parameter (energy units)\n";
    //    fo<<calculator.cfun.data1D[0].GetDiracAmplitude().real()<<"\n";
    //}
    
    fo<<"BathNumberOfOscillators:\n";
    fo<<"# Number of environment spectral densities\n";
    fo<<calculator.mfun.GetSize()<<"\n";
    
    fo<<"BathSpectralDensities:\n";
    fo<<"# List of spectral densities\n";
    for(int ind=0; ind<calculator.mfun.GetSize();ind++)
        fo<<calculator.spdf.data1D[ind].GetTitle()<<"\n";
    
    // system-bath coupling matrices
    
    fo<<"BathSystemCorrelatedCouplingMagnitudes:\n";
    for(int ind=0; ind<calculator.mfun.GetSize();ind++)
    {
        fo<<"# offdiagonal s-b coupling correlation triangular matrix\n";
        for(int ia = 0; ia<calculator.numT; ia++)
            for(int ib = 0; ib<=ia; ib++)
            {
                fo<< calculator.kij.data3D[ia][ib][ind];

                if(ib == ia) fo<<"\n";
                else fo<<"\t";
            }

        fo<<"# diagonal  s-b coupling correlation triangular matrix\n";
        for(int ia = 0; ia<calculator.numT; ia++)
            for(int ib = 0; ib<=ia; ib++)
            {
                fo<< calculator.gij.data3D[ia][ib][ind];

                if(ib == ia) fo<<"\n";
                else fo<<"\t";
            }
    }


    
    return 0;
}
