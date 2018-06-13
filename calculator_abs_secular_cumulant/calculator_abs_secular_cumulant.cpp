#include"../toolsFFT/toolsFFT.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../storage/storage.hpp"
#include"../constants/constants.hpp"
#include"../numericalSD/numericalSD.hpp"
#include"../feinman2sideddiagram1/feinman2sideddiagram1.hpp"
#include"../calculator_redfield/calculator_redfield.hpp"


#include"calculator_abs_secular_cumulant.hpp"


calculator_abs_secular_cumulant::calculator_abs_secular_cumulant():
communicator_3rd()
{
    eigensys = true;
    
    outputstring = "# QCFP calculator_abs_secular_cumulant class signal\n";
    
    liouville_pathway[0] = 0; // this one turns the pathway study on or off
    liouville_pathway[1] = 0; // this one sets the ground state
    liouville_pathway[2] = 1; // this one sets the excited state
    
    noComplexLifetimes=0;
    markovian = 0;
}



calculator_abs_secular_cumulant::calculator_abs_secular_cumulant(string& ifname):
communicator_3rd()
{
    eigensys = true;
    
    outputstring = "# QCFP calculator_abs_secular_cumulant class signal\n";
    
    liouville_pathway[0] = 0; // this one turns the pathway study on or off
    liouville_pathway[1] = 0; // this one sets the ground state
    liouville_pathway[2] = 1; // this one sets the excited state
    
    noComplexLifetimes=0;
    markovian = 0;
    
    toolsIO tio;
    // reading the whole input file
    tio.ReadWholeFile(ifname,inputfile,1);
    
    //std::cout<<inputfile.str();
    //std::cout<<"this is supposed to be the end of file\n";
    //exit(0);
    
    ifstream istr(ifname.c_str());
    if(!istr.is_open())
    {
        cout<<"Error: input file cannot be read\nquitting.\n";
        return;
    }
    tio.StreamSkipTrailers(istr);
    string tstr; getline(istr,tstr);
    ReadSystem(istr);
    ReadBath(istr);
    ReadExperimentLin(istr);
    
    storage<int> iri(2);
    // specific methods for this calculator
    if(tio.LookUpAndReadSquareMatrix<int>(
        "MethodSelectPathway:",
        "Reading the pathway selection\n", 
        1,2, iri,inputfile.str()))
    {
        liouville_pathway[0]=1;
        liouville_pathway[1]=iri.data2D[0][0];
        liouville_pathway[2]=iri.data2D[0][1];
    }
    
    std::size_t found;

    found = experiment_complexity.find("NoComplexLifetimes");
    if(found!=std::string::npos)
	noComplexLifetimes = 1;

    found = experiment_complexity.find("Markovian");
    if(found!=std::string::npos)
	markovian = 1;

    found = experiment_complexity.find("NonMarkovian");
    if(found!=std::string::npos)
	markovian = 0;

}



calculator_abs_secular_cumulant::~calculator_abs_secular_cumulant()
{
}
// calculator_abs_secular_cumulant::calculator_abs_secular_cumulant(
// storage<double>& ielevs,
// storage<dvector3d>& iedips,
// storage<complexv>& idephasings
// )
// {
// /*
//     eigensys = true;
//     outputstring = "# QCFP calculator_abs_levels class signal\n";
//     // acquire system parameters
//     evals = ielevs;
//     edips = iedips;
//     dephasings = idephasings;
//     
//     averagingtype = 0;
//     vecE1 = dvector3d(0,0,1);
//     vecE2 = dvector3d(0,0,1);
//     
//     liouville_pathway[0] = 0; // this one turns the pathway study on or off
//     liouville_pathway[1] = 0; // this one sets the ground state
//     liouville_pathway[2] = 1; // this one sets the excited state
//     noComplexLifetimes=0;
// */
// }
// 
// void calculator_abs_secular_cumulant::AddLineshapes
// (
//  storage<asymptoticLF_complexv>& igfun,
//  storage<asymptoticLF_complexv>&  imfun,
//  storage<double>& igij,
//  storage<double>& ikij,
//  double itempr)
// {
// /*
//     gfun = igfun;
//     mfun = imfun;
//     gij = igij; // amplitudes for all pairs
//     kij = ikij;
//     tempr = itempr;
// */
// }
// 


/*     
 * interpolationF<complexv> calculator_abs_secular_cumulant::Launch(double iifre, double iffre, int inump)
 * {
 *    if(ready == 0)
 *    {
 *        cout<<"Error: calculator_abs_levels::Launch cannot proceed\n";
 *    }
 *    
 *    ifre1 = iifre;
 *    ffre1 = iffre;
 *    nump = inump;
 *    
 *    Launch();
 *    
 *    // output
 *    return signal;
 * 
 * }
 * 
 */



void calculator_abs_secular_cumulant::Launch()
{
    //if(ready == 0)
    //{
    //    cout<<"Error: calculator_abs_secular_cumulant::Launch cannot proceed\n";
    //}
    cout<<"# Lauching calculator_abs_secular_cumulant::Launch\n";
    
    
    constants cst;
    double cfre = 0.5*(ifre1+ffre1);
    
    double freSt = (ffre1-ifre1)/nump;
    double freF = ffre1-ifre1;
    
    complexv* idat=0;
    if(markovian)
    {
        // will calculate direct frequency domain signal
        signal1d = interpolationF<complexv>(ifre1,freSt,nump);
        idat = signal1d.DirectAccessD();
        
    }
    
    
    // making time variables
    int timeN = nump/2;
    double timeS = cst.pi2/freF;
    double timeF = timeS*timeN;
    
    //making response functions (these will be used only in cumulant case
    interpolationF<complexv> asymrest(0,timeS,timeN);
    interpolationF<complexv> asymres(0,timeS,timeN);
    
    // system characteristics
    int& numg = numG;
    int& nume = numE;
    
    int numBosc =mfun.GetSize();
    
    //cout<<"somewhere?\n";
    
    // creating eigenstates (if missing)
    if(!evals.IsAlloc())
        MakeESSystem();
    FinalizeESSystem(!markovian);
    
/*    if(true)
    {
        calculator_redfield calcR(evals);
        calcR.AddFluctuations(mfun,kij,gij);
        
        calcR.AddTransportRates();
        // next I make dephasings and reorganizations
        if(!dephasings.IsAlloc())
        {
            dephasings.Allocate(numT,numT);
            dephasings = calcR.AddLifetimeDephasings();
            if(markovian)
                dephasings = calcR.AddPureDephasings();
        }
        if(!reorganizations.IsAlloc())
        {
            reorganizations.Allocate(numG+numE,numG+numE);
            reorganizations = calcR.GetReorganizations();
        }
        
 */       
        
        if(noComplexLifetimes)
        {
            int n2,n1;
            dephasings.GetSize(n2,n1);
            for(int j1=0; j1<n2; j1++)
                for(int j0=0; j0<n1;j0++)
                    dephasings.data2D[j1][j0] = dephasings.data2D[j1][j0].real();
                // notice that this overwrites dephasings.data2D as well
        }
        
    //}
    
    
    
    int iig = 0;
    int ifg = numg;
    int iie = numg;
    int ife = (nume+numg);
    if(liouville_pathway[0])
    {
        iig = liouville_pathway[1]; // this one sets the ground state
        ifg = iig+1;
        iie = liouville_pathway[2]; // this one sets the excited state
        ife = iie+1;
    }
    
    if(iig>numg || ifg>numg || iie > numg+nume || ife > numg+nume)
    {
        cout<<"Error: wrong Liouville pathway\n";
        return;
    }
    
    
    // lab configuration
    dvector3d ies[2];
    dvector3d iks[2];
    double ios[2];
    
    // polarizations
    ies[0] = vecE1;
    ies[1] = vecE2;
    
    // wavevectors : these are laser directions (not of interactions) : arbitrary units - they are normalized later anyway
    iks[0] = veck1;
    iks[1] = veck2;
    
    // thse will be fourier frequencies . those of interactions +k = +w; -k = -w.
    ios[0] = xomega1;
    ios[1] = xomega2;
    
    int beyonddipole =0;
    if(eten.IsAlloc() || ed_m.IsAlloc()) beyonddipole =1;
    interaction transition(2, beyonddipole,averagingtype);
    transition.PopulateSys(edips.data2D,eten.data2D,ed_m.data2D);
    ios[0] = xomega1;
    ios[1] = -xomega2;
    transition.PopulateFields(ies,iks,ios);
    
    
    
    // sum over energy levels
    for(int ig = iig; ig< ifg; ig++)
        for(int ie = iie; ie< ife; ie++)
        {
            transition.AssignLevels(0,ig,ie);
            transition.AssignLevels(1,ig,ie);
            
            
            if(markovian)
            {
                // direct calculation of the markovian Lorentzian
                complexv damplitude = transition.GetAveragedAmplitude();
                for(int ind1 = 0; ind1<nump; ind1++)
                {
                    double f1 = ifre1 + freSt*ind1;
                    complexv cval  = (f1 - (-2*(reorganizations.data2D[ie][ig]-reorganizations.data2D[ig][ig]) +evals.data1D[ie]-evals.data1D[ig]-coni*(dephasings.data2D[ie][ig]+naturallinewidth)));
                    cval =  grPops.data1D[ig]/cval;
                    cval *= -damplitude;
                    idat[ind1] += cval;
                }
            }
            else
            {
                // cumulant expansion case
                complexv omegas = -cfre + ( -2*(reorganizations.data2D[ie][ig]-reorganizations.data2D[ig][ig]) + evals.data1D[ie]-evals.data1D[ig]-coni*(dephasings.data2D[ie][ig]+naturallinewidth));
                feinman2sideddiagram1 diagram( omegas, transition);
                diagram.assignLineshape(gfun.data1D, gij.data3D[ie][ie],gij.data3D[ig][ig],gij.data3D[ig][ie],gij.data3D[ie][ig],numBosc);
                diagram.propagate(asymrest);
                asymrest *= (grPops.data1D[ig]);
                asymres += asymrest;
            }
        }
        
        if(markovian)
            return;
        
        // the rest is for time domain
        
        //cout<<"Error2\n";
        
        // shifting zero point due to Fourier transformation
        idat = asymres.DirectAccessD();
        idat[0] /= 2.0;
        idat = 0;
        //for(int ind = 0; ind< timeN; ind++)
        //{
        //    cout<<"sig: "<<idat[ind]<<"\n";
        //}
        
        //    cout<<asymres.Get(0)<<"\n";
        //    cout<<asymres.Get(0.1)<<"\n";
        
        
        // doing Fourier
        interpolationF<complexv> asymresFFTi(0,timeS,nump);
        interpolationF<complexv> asymresFFTf(0,(ffre1-ifre1)/nump,nump);
        //asymresFFTi *= complexv(0,0);
        asymresFFTi.DirectAssign(asymres);
        
        //    cout<<asymresFFTi.Get(0)<<"\n";
        //    cout<<asymresFFTi.Get(0.1)<<"\n";
        //for(int ind = 0; ind< nump; ind++)
        //{
        //    cout<<"sig: "<<asymresFFTi.Get(timeS*ind)<<"\n";
        //}
        
        
        toolsFFT fft;
        //asymresFFTf = asymresFFTi;
        asymresFFTf = fft.executeN(asymresFFTi);
        
        //for(int ind = 0; ind< nump; ind++)
        //{
        //    cout<<"sig: "<<(ffre+ifre)/2+(ffre-ifre)/nump*ind<<"\t"<<asymresFFTf.Get((ffre-ifre)/nump*ind)<<"\n";
        //}
        
        
        
        //    cout<<asymresFFTf.Get(0)<<"\n";
        //    cout<<asymresFFTf.Get(0.1)<<"\n";
        
        // shifting frequencies
        fft.SwapSides(asymresFFTf);
        asymresFFTf.UpdateAxis(ifre1,(ffre1-ifre1)/nump);
        
        //for(int ind = 0; ind< nump; ind++)
        //{
        //    cout<<"sig: "<<ifre+(ffre-ifre)/nump*ind<<"\t"<<asymresFFTf.Get(ifre+(ffre-ifre)/nump*ind)<<"\n";
        //}
        
        //    cout<<asymresFFTf.Get(0)<<"\n";
        //    cout<<asymresFFTf.Get(0.1)<<"\n";
        
        // output 
        signal1d = asymresFFTf;
}

