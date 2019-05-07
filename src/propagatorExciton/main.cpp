#include<iostream>

#include"propagatorExciton.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../communicator_3rd/communicator_3rd.hpp"

using std::cin;
using std::cout;
using std::ifstream;
using std::ofstream;


// this programs solves the Redfield equation for excitons
// it propagates the Green's function specified on the command-line
// as two indices

int main(int argc, const char * argv[])
{
    toolsIO tio;

    string ifilename(argv[1]);
    string ofilename(argv[2]);
    string str="";


    if (argc!=3)
    {
            cout<<"Error: specify the input file and the output file\n";
            return 0;
    }

    ifstream istr(ifilename.c_str());
    if(!istr.is_open())
    {
        cout<<"Error: input file cannot be read\n";
	return 0;
    }
    communicator_3rd reader;
	tio.ReadWholeFile(ifilename,reader.inputfile,1);

    tio.StreamSkipTrailers(istr);
    string tstr; getline(istr,tstr);
    reader.ReadSystem(istr);
    reader.ReadBath(istr);

    // next reading specific information for the propagations.

    string keyword;
    cout<<"Reading text-format keyword denoting approach of propagation:\n";
    // keyword denoting approach of propagation
    // available options are:
    // SecularRedfield
    // MarkovianRedfield
    // LindbladRedfield
    str = "";
    tio.StreamSkipTrailers(&istr);
    getline(istr,str);
    tio.StreamSkipTrailers(str);
    keyword = str;
    cout<<"Got: "<<str<<"\n";

    // excitonic block under consideration:
    // 10, 20, 11, 21, 22
	storage<double> arr(2);
	storage<int> iri(2);
    int iblock;
	// specific methods for this calculator
	if(tio.LookUpAndReadSquareMatrix<int>(
			"MethodBlockSelector:",
			"Reading the block to propagate\n",
			1,1, iri,reader.inputfile.str()))
		{
			iblock=iri.data2D[0][0];
		}

    // setting up number of elements
    int numL, numR;
    // setting up number of elements
    if(iblock==0)
    {
      numL = 1;
      numR = 1;
	//shL = 0;
	//shR = 0;
      reader.band = 0;
    }
    else if(iblock==10)
    {
      numL = reader.numE;
      numR = 1;
	//shL = 1;
	//shR = 0;
      reader.band = 1;
    }
    else if(iblock==20)
    {
      numL = reader.numF;
      numR = 1;
	//shL = 1+reader.numE;
	//shR = 0;
      reader.band = 2;
    }
    else if(iblock==11)
    {
      numL = reader.numE;
      numR = reader.numE;
	//shL = 1;
	//shR = 1;
      reader.band = 1;
    }
    else if(iblock==21)
    {
      numL = reader.numF;
      numR = reader.numE;
	//shL = 1+reader.numE;
	//shR = 1;
      reader.band = 2;
    }
    else if(iblock==22)
    {
      numL = reader.numF;
      numR = reader.numE;
	//shL = 1+reader.numE;
	//shR = 1+reader.numE;
      reader.band = 2;
    }
    else
    {
        cout<<"Error: exciton block specified incorrectly\n";
        return 0;
    }



    // Initial condition matrix: real parts
    storage<complexv> initials(2);
    initials.Allocate(numL,numR);
	if(tio.LookUpAndReadSquareMatrix<double>(
			"MethodInitialConditionMatrixRe:",
			"Reading the initial condition matrices: Re\n",
			numL,numR, arr,reader.inputfile.str()))
		{
            for(int ind=0; ind<numL; ind ++)
            for(int inz=0; inz<numR; inz ++)
            {
                initials.data2D[ind][inz]=arr.data2D[ind][inz];
            }
        }
	if(tio.LookUpAndReadSquareMatrix<double>(
			"MethodInitialConditionMatrixIm:",
			"Reading the initial condition matrices: Im\n",
			numL,numR, arr,reader.inputfile.str()))
		{
            for(int ind=0; ind<numL; ind ++)
            for(int inz=0; inz<numR; inz ++)
            {
                initials.data2D[ind][inz]=complexv(initials.data2D[ind][inz].real(),arr.data2D[ind][inz]);
            }
        }

  double timestep=0.1;
	if(tio.LookUpAndReadSquareMatrix<double>(
			"MethodOutputTimeStep:",
			"Reading the output timestep\n",
			1,1, arr,reader.inputfile.str()))
		{
            timestep = arr.data2D[0][0];
        }

  int timepoints=0;
	if(tio.LookUpAndReadSquareMatrix<int>(
			"MethodOutputTimePoints:",
			"Reading the number of points on the output\n",
			1,1, iri,reader.inputfile.str()))
		{
            timepoints = iri.data2D[0][0];
        }


	//implementation of internaltimeS and internaltimeN
  double internaltimestep=0;
	if(tio.LookUpAndReadSquareMatrix<double>(
    			"MethodInternalTimeStep:",
    			"Reading the timestep for correlation function\n",
    			1,1, arr,reader.inputfile.str()))
  		{
          internaltimestep = arr.data2D[0][0];
     	 }

  int internaltimepoints=0;
	if(tio.LookUpAndReadSquareMatrix<int>(
    			"MethodInternalTimePoints:",
    			"Reading the number of points of correlation function\n",
    			1,1, iri,reader.inputfile.str()))
  		{
          internaltimepoints = iri.data2D[0][0];
      	}


////////////////////////////////////////////////////////////////////

    std::size_t found;
    //found = keyword.find("NonSecular");
    //if(found!=std::string::npos)
	//reader.nonsecular = 1;
    // secular holds only for Markovian
    // nonsecular simple = full Redfield

    reader.MakeESSystem();
    propagatorExciton pE(reader);

    found = keyword.find("NonSecular");
    if(found!=std::string::npos){
	pE.flagNonsecular=1;
    }

    found = keyword.find("Lindblad");
    if(found!=std::string::npos){
	pE.flagMarkovian=1;
	pE.flagNonsecular=1;
	pE.calcR->flagLindblad = 1;
    //pE.calcR->nonsecular=1
    }


    found = keyword.find("ModRed");
    if(found!=std::string::npos){
	pE.flagModred = 1;
	pE.flagMarkovian=1;
	pE.flagNonsecular=0;
	pE.calcR->flagLindblad = 0;
    //pE.calcR->nonsecular=0
    }

		found = keyword.find("NonMarkovian");
    if(found!=std::string::npos){
			pE.flagMarkovian=0;
			pE.flagMemoryWithCfun=0;
			pE.flagNonsecular=1;
    }
		found = keyword.find("NonMarkovianCFun");
		if(found!=std::string::npos){
			pE.flagMarkovian=0;
			pE.flagMemoryWithCfun=1;
			pE.flagNonsecular=1;
    }
		found = keyword.find("TC2");
    if(found!=std::string::npos){
			pE.flagMemoryWithCfun = 0+pE.flagMemoryWithCfun;
    }
		found = keyword.find("TCL2");
    if(found!=std::string::npos){
			pE.flagMemoryWithCfun = 10+pE.flagMemoryWithCfun;
    }

     // setting initial condition
    pE.dmatrix0 = initials;
    pE.internaltimeS = internaltimestep;
    pE.internaltimeN =internaltimepoints;

    // making time trace variable
    storage<double> times(1);
    times.FillLinear(0,timestep,timepoints);

    storage<complexv> dmt(3);
    std::cout<<"Calculating\n";

    // checking the block for propagation
    pE.SetBlock(iblock);
		pE.manifold0End = 0;
    dmt = pE.PropagateDM(times);




    // printing result to output file:
    ofstream ofs(ofilename.c_str());

    // final output
		ofs.precision(16);
    ofs<<"# time"<<"\t";
        for(int ia = 0; ia< numL; ia++)
        for(int ib = 0; ib< numR; ib++)
        {
            ofs<<"("<<ia<<" "<<ib<<")re\t";
            ofs<<"("<<ia<<" "<<ib<<")im";
            if(ia*numR+ib == numL*numR-1) ofs<<"\n";
            else ofs<<"\t";
        }
        for(int it = 0; it< timepoints; it++)
        {
            ofs<<it*timestep<<"\t";
            for(int ia = 0; ia< numL; ia++)
            for(int ib = 0; ib< numR; ib++)
            {
                ofs<<dmt.data3D[it][ia][ib].real()<<"\t";
                ofs<<dmt.data3D[it][ia][ib].imag();
                if(ia*numR+ib == numL*numR-1) ofs<<"\n";
                else ofs<<"\t";
            }
        }




}
