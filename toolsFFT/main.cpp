#include"toolsFFT.hpp"

#include"../storage/storage.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../interpolationF/interpolationF.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"

#include"../toolsIO/toolsIO.hpp"

#include<iostream>
#include<fstream>
#include<string>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::string;


int main(int argc, char **argv)
{
    
	if(argc == 1)
	{
		cout<<"This is a small program to calculate Fourier transform\n";
		cout<<"for a complex function (two-column dataset).\n";
		cout<<"\n";
		cout<<"    Usage is controlled by the command-line arguments:\n";
		cout<<"1-st argument is the file-name to read\n";
		cout<<"2-nd argument is the number of values to read\n";
		cout<<"3-rd argument is the column number to read as real\n";
		cout<<"4-th argument is the column number to read as imaginary.\n";
		cout<<"\n";
		return 0;
	}
	
	if(argc !=5)
	{
		cout<<"Wrong number of Arguments\n";
		return 0;
	}


	// proceeding
    toolsIO tio;


	string* sargv = new string[argc];
	for(int ind = 0; ind<argc; ind++) 
		sargv[ind] = string(argv[ind]);


	// number of points
	int Num = tio.fromString<int>(sargv[2]);
	int cRe = tio.fromString<int>(sargv[3]) - 1; // in c++ counting from zero
	int cIm = tio.fromString<int>(sargv[4]) - 1;
	// all reading is done


	// reading dataset from file
	ifstream ifs(sargv[1].c_str());
	if(!ifs.is_open())
	{
		cout<<"Error: toolsFFT: cannot open file!\n";
		return 0;
	}

    
    // preparing buffers for source data
    complexv* dinp = new complexv[Num];
    complexv* dout = new complexv[Num];
    
    
	string singleLine;
    int counter = 0;

    while(!ifs.eof())
    {
        tio.StreamSkipTrailers(&ifs);
        getline(ifs, singleLine);
        
        // beaking line into words
        int numw;
        string* words;
        words = tio.toWords(numw,singleLine);
        
        // reading numbers
        if(numw >cRe && numw> cIm && counter <Num)
        {
            dinp[counter]=complexv( tio.fromString<double>(words[cRe]),tio.fromString<double>(words[cIm]) );
            counter++;
        }
        
        delete[] words;
    }
    // closing file
    ifs.close();

    // data ready and making FFT:
    toolsFFT tff;
    tff.executeP(dout,dinp,Num);

    // FFT done
    // outputing data
    ofstream ofs((sargv[1]+".fft").c_str());    
    for(int ind = 0; ind< Num; ind++)
    {
        ofs<<dout[ind].real()<<"\t"<<dout[ind].imag()<<"\n";
    }
    ofs.close();
    
    // deleting datasets
    delete[] sargv;
    delete[] dinp;
    delete[] dout;

}
