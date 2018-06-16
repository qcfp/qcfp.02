#pragma once
#include"../complexv/complexv.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"

class correlation4points
{
private:

    // no allocation is performed in this object: all parameters are copies
    // no deallocation is done
    
        // bath:
        asymptoticLF_complexv* gf;
        //double** ampl;
        int sa;
        int sb;
        int sc;
        int sd;
        int se;
        int se1;
        
        int ready;
        
        
        // for the multi-mode bath
        int numO;
        double*** amplt;

public:
        correlation4points();
        correlation4points(asymptoticLF_complexv* gij,double*** gijAmpl,int inumO, int id,int ic, int ib, int ia);
        //correlation4points(asymptoticLF_complexv* gij,double** gijAmpl, int id,int ic, int ib, int ia);
        ~correlation4points();


complexv Get(double& t4,double& t3,double& t2,double& t1);

complexv GetIncoherent(double& t3,double& t2,double& t1);


};

