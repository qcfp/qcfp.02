#include"correlation4points.hpp"

correlation4points::correlation4points()
{
    // do not much
    gf = 0;
    amplt = 0;
    
    ready = 0;
}

correlation4points::~correlation4points()
{
    // do not much
    gf = 0;
    amplt = 0;
    
    ready = 0;
}


//correlation4points::correlation4points(asymptoticLF_complexv* gij,double** gijAmpl,int isd,int isc, int isb, int isa)
//{
//    // this is for the single bath mode
//    
//    // setting up phase for coherent diagram
//    numO = 1;
//    gf = gij;
//    ampl = gijAmpl;    
//    ready = 1;
//    sa=isa;
//    sb=isb;
//    sc=isc;
//    sd=isd;
//    
//    // the following are to be used with the incoherent secular transport
//    se = sa;
//    se1 = sd;
//}

correlation4points::correlation4points(asymptoticLF_complexv* gij,double*** gijAmpl,int inumO, int isd,int isc, int isb, int isa)
{
    // this is for the multi-mode bath mode
    
    // setting up phase for coherent diagram
    gf = gij;
    amplt = gijAmpl;
    numO = inumO;
    ready = 1;
    sa=isa;
    sb=isb;
    sc=isc;
    sd=isd;
    
    // the following are to be used with the incoherent secular transport
    se = sa;
    se1 = sd;
}


complexv correlation4points::Get(double& t4, double& t3, double& t2, double& t1)
{
    
    //cout<<"t4 = "<<t4<<"\n";
    //cout<<"t3 = "<<t3<<"\n";
    //cout<<"t2 = "<<t2<<"\n";
    //cout<<"t1 = "<<t1<<"\n";
    
    //cout<<"dcba = "<<sd<<sc<<sb<<sa<<"\n";
    
    if(ready == 0)
    {
        cout<<"Error in correlation4points\n";
        return 0;
    }
    
    complexv fval;
    complexv cvalue(0.0,0.0);
    
     //       complexv ggg = complexv(0.1,-1.0);


    for(int io = 0; io< numO; io++)
    {
        //if(numO > 1)
        //{
        //        gf = gft[io];
                //ampl = amplt[io];
        //}
    
        double ampl_sd__sd_ = amplt[sd][sd][io];
        double ampl_sd__sc_ = amplt[sd][sc][io];
        double ampl_sa__sd_ = amplt[sa][sd][io];
        double ampl_sa__sc_ = amplt[sa][sc][io];
        double ampl_sc__sc_ = amplt[sc][sc][io];
        double ampl_sd__sb_ = amplt[sd][sb][io];
        double ampl_sc__sb_ = amplt[sc][sb][io];
        double ampl_sb__sb_ = amplt[sb][sb][io];
        double ampl_sc__sa_ = amplt[sc][sa][io];
        double ampl_sb__sa_ = amplt[sb][sa][io];
        double ampl_sa__sa_ = amplt[sa][sa][io];
        double ampl_sb__sd_ = amplt[sb][sd][io];
        double ampl_sa__sb_ = amplt[sa][sb][io];
        double ampl_sd__sa_ = amplt[sd][sa][io];
        
        
//        cout<<"A dd = "<<ampl_sd__sd_<<"\n";
//        cout<<"A dc = "<<ampl_sd__sc_<<"\n";
//        cout<<"A ad = "<<ampl_sa__sd_<<"\n";
//        cout<<"A ac = "<<ampl_sa__sc_<<"\n";
//        cout<<"A cc = "<<ampl_sc__sc_<<"\n";
//        cout<<"A db = "<<ampl_sd__sb_<<"\n";
//        cout<<"A cb = "<<ampl_sc__sb_<<"\n";
//        cout<<"A bb = "<<ampl_sb__sb_<<"\n";
//        cout<<"A ca = "<<ampl_sc__sa_<<"\n";
//        cout<<"A ba = "<<ampl_sb__sa_<<"\n";
//        cout<<"A aa = "<<ampl_sa__sa_<<"\n";
//        cout<<"A bd = "<<ampl_sb__sd_<<"\n";
//        cout<<"A ab = "<<ampl_sa__sb_<<"\n";
//        cout<<"A da = "<<ampl_sd__sa_<<"\n";

        
        fval = gf[io].Get(t4-t3);
        //fval = (t4>t3) ? ggg*(t4-t3) : ggg.conjugate()*(t3-t4);
    //cout<<"gfun: "<<io<<" "<<t4-t3<<" "<<fval<<"\n";
        cvalue += (-ampl_sd__sd_+ampl_sd__sc_)*fval+conj((+ampl_sa__sd_-ampl_sa__sc_)*fval);

        fval = gf[io].Get(t3-t2);
        //fval = (t3>t2) ? ggg*(t3-t2) : ggg.conjugate()*(t2-t3);
    //cout<<"gfun: "<<io<<" "<<t3-t2<<" "<<fval<<"\n";
       cvalue += (-ampl_sc__sc_+ampl_sd__sc_-ampl_sd__sb_+ampl_sc__sb_)*fval;
        
        fval = gf[io].Get(t2-t1);
        //fval = (t2>t1) ? ggg*(t2-t1) : ggg.conjugate()*(t1-t2);
    //cout<<"gfun: "<<io<<" "<<t2-t1<<" "<<fval<<"\n";
        cvalue += (-ampl_sb__sb_+ampl_sc__sb_-ampl_sc__sa_+ampl_sb__sa_)*fval;

        fval = gf[io].Get(t1-t4);
        //fval = (t1>t4) ? ggg*(t1-t4) : ggg.conjugate()*(t4-t1);
    //cout<<"gfun: "<<io<<" "<<t1-t4<<" "<<fval<<"\n";
        cvalue += (-ampl_sa__sa_+ampl_sb__sa_)*fval+conj((-ampl_sb__sd_+ampl_sa__sd_)*fval);
        
        fval = gf[io].Get(t4-t2);
        //fval = (t4>t2) ? ggg*(t4-t2) : ggg.conjugate()*(t2-t4);
    //cout<<"gfun: "<<io<<" "<<t4-t2<<" "<<fval<<"\n";
        cvalue += (-ampl_sd__sc_+ampl_sd__sb_)*fval+conj((+ampl_sa__sc_-ampl_sa__sb_)*fval);

        fval = gf[io].Get(t3-t1);
        //fval = (t3>t1) ? ggg*(t3-t1) : ggg.conjugate()*(t1-t3);
    //cout<<"gfun: "<<io<<" "<<t3-t1<<" "<<fval<<"\n";
        cvalue += (+ampl_sd__sb_-ampl_sd__sa_-ampl_sc__sb_+ampl_sc__sa_)*fval;
    }

//      fval =  -100.0*ggg.conjugate()*(t3-t2);
//      fval += -100.0*ggg.conjugate()*(t1-t4);
//      cvalue = fval;
      

//    for(int io = 0; io< numO; io++)
//    for(int i2 = 0; i2<4; i2++)
//    for(int i1 = 0; i1<4; i1++)
//    cout<<io<<" "<<i2<<" "<<i1<<" "<<amplt[io][i2][i1]<<"\n";
        return cvalue;
}


complexv correlation4points::GetIncoherent(double& t3, double& t2, double& t1)
{
    if(ready == 0)
    {
        cout<<"Error in correlation4points\n";
        return 0;
    }

    complexv fval;
    complexv cvalue(0.0,0.0);

    for(int io = 0; io< numO; io++)
    {
        //if(numO > 1)
        //{
//                gf = gft[io];
                //ampl = amplt[io];
        //}
        complexv cval2(0.0,0.0);
    
        double ampl_se__se_ = amplt[se][se][io];
        double ampl_sb__sb_ = amplt[sb][sb][io];
        double ampl_sc__sb_ = amplt[sc][sb][io];
        double ampl_sc__sc_ = amplt[sc][sc][io];
        double ampl_sc__se1 = amplt[sc][se1][io];
        double ampl_sb__se1 = amplt[sb][se1][io];
        double ampl_sc__se_ = amplt[sc][se][io];
        double ampl_sb__se_ = amplt[sb][se][io];


//        cout<<io<<" "<<t3<<" "<<t2<<" "<<t1<<" ";//<<fval<<" ";

        
        fval = gf[io].Get(t1);
//        cout<<fval<<" ";
        
        cvalue += -ampl_se__se_*fval;

        fval = gf[io].Get(t3);
//        cout<<fval<<" ";

        cvalue += (-ampl_sb__sb_+ampl_sc__sb_)*fval+conj((-ampl_sc__sc_+ampl_sc__sb_)*fval);
        cval2 += +(-ampl_sc__se1+ampl_sb__se1)*fval;

        fval = gf[io].Get(t2);
//        cout<<fval<<" ";

        cvalue += (+ampl_sc__se_-ampl_sb__se_)*fval;
        cval2 += +(-ampl_sc__se1+ampl_sb__se1)*fval;

        fval = gf[io].Get(t1+t2);
//        cout<<fval<<" ";

        cvalue += (+ampl_sb__se_-ampl_sc__se_)*fval;

        fval = gf[io].Get(t2+t3);
//        cout<<fval<<" ";

        cvalue += (+ampl_sb__se_-ampl_sc__se_)*fval;
        cval2 += +(+ampl_sc__se1-ampl_sb__se1)*fval;

        fval = gf[io].Get(t1+t2+t3);
//        cout<<fval<<" ";

        cvalue += (-ampl_sb__se_+ampl_sc__se_)*fval;
        
        
        cvalue += complexv(0.0, 2*cval2.imag() );
        
//        cout<<cvalue<<"\n";
    }

        return cvalue;
}

