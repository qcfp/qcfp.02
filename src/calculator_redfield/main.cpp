#include<iostream>

#include"calculator_redfield.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../communicator_3rd/communicator_3rd.hpp"
#include"../numericalSD/numericalSD.hpp"

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
    int numL, numR, shL, shR, iblock;

    if (argc!=3)
    {
            cout<<"Error: specify the input file and the output file\n";
            return 0;
    }

    string ifilename(argv[1]);
    string ofilename(argv[2]);
    string str="";

    // perform calculations based on filename:
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
    cout<<"Reading text-format keyword denoting the type of superoperator:\n";
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

	storage<double> arr(2);
	storage<int> iri(2);
	// specific methods for this calculator
	if(tio.LookUpAndReadSquareMatrix<int>(
			"MethodBlockSelector:",
			"Reading the block to propagate\n", 
			1,1, iri,reader.inputfile.str()))
		{
			iblock=iri.data2D[0][0];
		}

    // block under consideration:
    // 10, 20, 11, 21, 22

    // setting up number of elements
    if(iblock==0)
    {
      numL = 1;
      numR = 1;
	shL = 0;
	shR = 0;
      reader.band = 0;
    }
    else if(iblock==10)
    {
      numL = reader.numE;
      numR = 1;
	shL = 1;
	shR = 0;
      reader.band = 1;
    }
    else if(iblock==20)
    {
      numL = reader.numF;
      numR = 1;
	shL = 1+reader.numE;
	shR = 0;
      reader.band = 2;
    }
    else if(iblock==11)
    {
      numL = reader.numE;
      numR = reader.numE;
	shL = 1;
	shR = 1;
      reader.band = 1;
    }
    else if(iblock==21)
    {
      numL = reader.numF;
      numR = reader.numE;
	shL = 1+reader.numE;
	shR = 1;
      reader.band = 2;
    }
    else if(iblock==22)
    {
      numL = reader.numF;
      numR = reader.numE;
	shL = 1+reader.numE;
	shR = 1+reader.numE;
      reader.band = 2;
    }
    else
    {
        cout<<"Error: exciton block specified incorrectly\n";
	return 0;
    }
        

////////////////////////////////////////////////////////////////////

    std::size_t found;
//    found = reader.experiment_complexity.find("Secular");
//    if(found!=std::string::npos)
//	reader.nonsecular = 0;
//    found = reader.experiment_complexity.find("NonSecular");
//    if(found!=std::string::npos)
//	reader.nonsecular = 1;

    reader.MakeESSystem();

    
   // SecularRedfield
    // MarkovianRedfield
    // LindbladRedfield
    //if(keyword.compare("MarkovianRedfield")==0)
    //pE.SetMarkovian(true);

	int flagMarkovian = 1;
    int flagNonsecular = 0;
    int flagLindblad = 0;
    int flagModred = 0;
    
    //found = reader.experiment_complexity.find("Markovian");
    //if(found!=std::string::npos)
	//flagMarkovian = 1;

    found = keyword.find("NonMarkovian");
    if(found!=std::string::npos){
	flagMarkovian = 0;
    flagNonsecular=1;
    }

    found = keyword.find("NonSecular");
    if(found!=std::string::npos){
	flagNonsecular=1;
    }

    found = keyword.find("Lindblad");
    if(found!=std::string::npos){
	flagMarkovian=1;
	flagNonsecular=1;
	flagLindblad = 1;
    }

    found = keyword.find("ModRed");
    if(found!=std::string::npos){
	flagModred = 1;
	flagMarkovian=1;
	flagNonsecular=0;
	flagLindblad = 0;
    }

    // specific part to the modified Redfield
    storage<asymptoticLF_complexv> g1un(1); // gfunctions first derivatives
    if(true){
    int numBosc = reader.mfun.GetSize();
    g1un.Allocate(numBosc);
    for(int ind=0; ind<numBosc; ind++)
    {
        numericalSD spd(reader.tempr,reader.spdf.data1D[ind]);
        g1un.data1D[ind] = spd.GetGfD1f();
    }
    }
    calculator_redfield calcR(reader.evals, reader.numG, reader.numE, reader.numF);
    calcR.AddFluctuations(reader.mfun,reader.cfun,reader.gfun,g1un,reader.kij,reader.gij,reader.zijkl);
        //calcR.AddFluctuations(reader.mfun,reader.cfun,reader.kij,reader.gij);
    calcR.tempr=reader.tempr;
	//calcR.nonsecular = reader.nonsecular;
    calcR.flagLindblad = flagLindblad;


    // Getting all required superoperators



	int num5,num4,num3,num2,num1;
	double deltaT=0;

    storage<complexv> superoperatorSG(2); // secular markovian dephasings
    storage<double> superoperatorSP(2); // secular markovian populations
    storage<complexv> superoperatorM(4); // markovian full
    storage<complexv> superoperatorR(5); // non-markovian full with memory

	storage<double> energiesReorg(2);
	energiesReorg =  calcR.GetReorganizations();

	if(flagNonsecular)
        calcR.AddMijkl(reader.zijkl);

	if(flagMarkovian)
	{
		if(flagNonsecular)
	        superoperatorM = calcR.GetRelaxationSuperoperatorM(iblock);
		else
		{
            if(flagModred)
			    superoperatorSP = calcR.GetTransportRatesModRed();
            else
			    superoperatorSP = calcR.AddTransportRates();
			calcR.AddLifetimeDephasings();
			superoperatorSG = calcR.AddPureDephasings();
		}
	}
    else // this is always nonsecular
	{
		num5 = 0;
       	superoperatorR = calcR.GetMemoryKernel(0, deltaT, num5 , iblock);
		superoperatorR.GetSize(num5,num4,num3,num2,num1);
	}



	

    
    // printing result to output file:
    ofstream ofs(ofilename.c_str());

	ofs.precision(12);

	ofs<<"# System eigenvector matrix: 1-exc eigenvectors in columns.\n";
    // print eigenvectors
    for(int indl=0; indl<reader.numE; indl++)
    {
        for(int inde=0; inde<reader.numE; inde++)
        {
            ofs<<reader.evec1.data2D[inde][indl];
            if(inde == reader.numE-1) ofs<<"\n";
            else  ofs<<"\t";
        }
    }


    // print eigenvectors F
    if(reader.band == 2) ofs<<"# 2-exc SSF (eigenvectors in columns):\n";
    if(reader.band == 2) 
    for(int indl=0; indl<reader.numE; indl++)
    {
        for(int indr=0; indr<=indl; indr++)
        {
            ofs<<indl<<"-"<<indr<<"\t";
            for(int inde=0; inde<reader.numF; inde++)
            {
                cout<<reader.evec2.data3D[inde][indl][indr];
                if(inde == reader.numF-1) ofs<<"\n";
                else  ofs<<"\t";
            }

        }
    }


	if(flagNonsecular)
	{

		if(!flagMarkovian)
		{
			// here is full Redfield relaxation operator
			ofs<<"# Eigenstate basis. Full nonsecular Markovian Redfield relaxation operator. Interaction picture.\n";

			ofs<<"# First energy gaps:\n";
			for(int id2=0;id2<numL;id2++)
			for(int id1=0;id1<numR;id1++)
				ofs<<id2<<" "<<id1<<" "<<reader.evals.data1D[id2+shL]-reader.evals.data1D[id1+shR]<<"\n";


			superoperatorR.GetSize(num5,num4,num3,num2,num1);

			for(int it=0; it<num5; it++)
			for(int id4=0;id4<num4;id4++)
			for(int id3=0;id3<num3;id3++)
			for(int id2=0;id2<num2;id2++)
			for(int id1=0;id1<num1;id1++)
				ofs<<-it*deltaT<<" "<<id4<<" "<<id3<<" "<<id2<<" "<<id1<<" "<<superoperatorR.data5D[it][id4][id3][id2][id1]<<"\n";

		}
		else
		{
			ofs<<"# Eigenstate basis. Full nonsecular Markovian Redfield relaxation operator.\n";

            superoperatorM = calcR.GetRelaxationSuperoperatorM(iblock);


	        for(int il2 = 0; il2<numL;il2++)
        	for(int il1 = 0; il1<numR;il1++)
        	{
        	    int il = il2*numR+il1;
	            for(int ir2 = 0; ir2<numL;ir2++)
	            for(int ir1 = 0; ir1<numR;ir1++)
	            {

			complexv cvalue(0,0);

//	                if(il2 == ir2 && il1 == ir1){
//	                    cvalue +=  cnni*(reader.evals.data1D[il2+shL]-energiesReorg.data2D[il2+shL][il2+shL]);
//	                    cvalue +=  coni*(reader.evals.data1D[il1+shR]-energiesReorg.data2D[il1+shR][il1+shR]);
//	                }
	                cvalue += superoperatorM.data4D[il2][il1][ir2][ir1];

			ofs<<il2<<" "<<il1<<" "<<ir2<<" "<<ir1<<" "<<cvalue.real()<<" "<<cvalue.imag()<<"\n";
	            }
	        }

		}


	}
	else
	{
		// secular markovian
		ofs<<"# Eigenstate basis. Secular Markovian Redfield relaxation operator.\n";
		

		ofs<<"# Population block:\n";
		for(int il=0;il<numL;il++)
		for(int ir=0;ir<numR;ir++)
		{
			ofs<<superoperatorSP.data2D[il+shL][ir+shL];
			if(ir == numR-1) ofs<<"\n"; else ofs<<"\t";
		}
    
		ofs<<"# Dephasing block (real parts):\n";
		for(int il=0;il<numL;il++)
		for(int ir=0;ir<numR;ir++)
		{
			ofs<<superoperatorSG.data2D[il+shL][ir+shL].real();
			if(ir == numR-1) ofs<<"\n"; else ofs<<"\t";
		}

		ofs<<"# Dephasing block (imaginary parts):\n";
		for(int il=0;il<numL;il++)
		for(int ir=0;ir<numR;ir++)
		{
			ofs<<superoperatorSG.data2D[il+shL][ir+shL].imag();
			if(ir == numR-1) ofs<<"\n"; else ofs<<"\t";
		}
	}
}
