#include<iostream>

#include"propagatorExcitonWF.hpp"
#include"../toolsIO/toolsIO.hpp"
#include"../communicator_3rd/communicator_3rd.hpp"

#include<iostream>


using std::cin;
using std::cout;
using std::ifstream;
using std::ofstream;


// this programs propagates the Davydov wavefunction for excitons
// it returns the electronic population as a function of time starting from
// a given electronic initial condition 
// and a thermal nuclear initial state
// of a given temperature

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

    // perform calculations based on filename:
    communicator_3rd reader;
	tio.ReadWholeFile(ifilename,reader.inputfile,1);
    
     tio.StreamSkipTrailers(istr);
    string tstr; getline(istr,tstr);
    reader.ReadSystem(istr);
    reader.ReadBath(istr);

    storage<double> arr(2);
    storage<int> iri(2);
    // next reading specific information for the propagations.

    string keyword;
    cout<<"Reading text-format keyword denoting approach of propagation:\n";
    // keyword denoting approach of propagation
    // available options are: 

    str = "";
    tio.StreamSkipTrailers(&istr);
    getline(istr,str);
    tio.StreamSkipTrailers(str);
    keyword = str;
    cout<<"Got: "<<str<<"\n";

    // excitonic block under consideration is 1:

 
    // Initial condition vector: real parts 1-st column, imaginary parts 2-nd column
    storage<complexv> initials(1);
    initials.Allocate(reader.numE);
	if(tio.LookUpAndReadSquareMatrix<double>(
			"MethodInitialConditionVector:",
			"Reading the initial condition vector: Re / Im\n", 
			reader.numE,2, arr,reader.inputfile.str()))
		{
            for(int ind=0; ind<reader.numE; ind ++)
            {
                initials.data1D[ind]=complexv(arr.data2D[ind][0],arr.data2D[ind][1]);
            }
        }

    double timestep=0.1;
	if(tio.LookUpAndReadSquareMatrix<double>(
			"MethodTimestep:",
			"Reading the output timestep\n", 
			1,1, arr,reader.inputfile.str()))
		{
            timestep = arr.data2D[0][0];
        }

    int timepoints;
	if(tio.LookUpAndReadSquareMatrix<int>(
			"MethodTimePoints:",
			"Reading the number of points on the output\n", 
			1,1, iri,reader.inputfile.str()))
		{
            timepoints = iri.data2D[0][0];
        }

////////////////////////////////////////////////////////////////////
	reader.MakeESSystem();
	wavefunctionD1 objd1(reader);
	propagatorExcitonWF pE(timestep);

     // setting initial condition
	objd1.MakeThermal();
	objd1.SetAlpha1(initials);
    	
    	pE.SwitchManifold(1,objd1);
    std::cout<<"Calculating\n";
    
    // printing result to output file:
    ofstream ofs(ofilename.c_str());

    // propagation and final output
	storage<double> times(1);
	times.FillLinear(0,timestep,timepoints);

        for(int it = 0; it< timepoints; it++)
        {
		double time = times.data1D[it];

		pE.PropagateWF(time,objd1);

            ofs<<time<<"\t";

        for(int ia = 0; ia< reader.numE; ia++)
        {
                ofs<<(objd1.alpha1[ia]).norm();
                if(ia == reader.numE-1) ofs<<"\n";
                else ofs<<"\t";
        }
	}


    
}
