#include<iostream>

#include"interpolationF.hpp"

#include"../complexv/complexv.hpp"
#include"../storage/storage.hpp"
#include"../toolsIO/toolsIO.hpp"

using std::cin;
using std::cout;


#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

using std::string;
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;

#include"../storage/storage.hpp"
#include"../interpolationF/interpolationF.hpp"

int main(int argc, const char * argv[])
{
	// this is binary to text and
	// vice versa conversion for complex and double numbers

	// usage:
	// program_name key input_file output_file
    // keys: t2b-double b2t-double t2b-complex b2t-complex
	// where t2b stands text-to-double, etc.
	/////////////////

//----------------------------------------------------
    if(argc !=4)
    {


        cout<<"Error: use three parameters: key ifile ofile\n";

        cout<<"keys: t2b-double b2t-double t2b-complex b2t-complex\n";
        exit(0);
    }

    string key(argv[1]);
    string fi1(argv[2]);
    string fi2(argv[3]);

    
    if(key == "t2b-double")
    {
       interpolationF<double> funct;
       funct.ReadF(fi1);
       funct.SaveF(fi2);
    }
    else if(key == "b2t-double")
    {
    	interpolationF<double> funct;
        funct.ReadF(fi1);
        funct.SaveFtxt(fi2);
    }
    else if(key == "t2b-complex")
    {
    	interpolationF<complexv> funct;
        funct.ReadF(fi1);
        funct.SaveF(fi2);
    }
    else if(key == "b2t-complex")
    {
    	interpolationF<complexv> funct;
        funct.ReadF(fi1);
        funct.SaveFtxt(fi2);
    }
    else
    {

        cout<<"Error: specified key is not recognized\n";

        cout<<"use: t2b-double b2t-double t2b-complex b2t-complex\n";
    }
    
     return 0;
}


