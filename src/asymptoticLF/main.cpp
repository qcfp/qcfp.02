
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
#include"../dvector3d/dvector3d.hpp"
#include"../interpolationF/interpolationF.hpp"

#include"asymptoticLF.hpp"


int main(int argc, const char* argv[])
{

//----------------------------------------------------
    if(argc !=4)
    {
        
        
        cout<<"Error: use three parameters: key ifile ofile\n";
        
        cout<<"keys: \n";
         cout<<"t2b-double\n";
         cout<<"b2t-double\n";
         cout<<"t2b-complex\n";
         cout<<"b2t-complex\n";
         cout<<"t2b-complexdouble\n";
         cout<<"b2t-complexdouble\n";
        exit(0);
    }

    string key(argv[1]);
    string fi1(argv[2]);
    string fi2(argv[3]);

    
    if(key == "t2b-double")
    {
       asymptoticLF<double> funct;
       funct.ReadF(fi1);
       funct.SaveF(fi2);
    }
    else if(key == "b2t-double")
    {
        asymptoticLF<double> funct;
        funct.ReadF(fi1);
        funct.SaveFtxt(fi2);
    }
    else if(key == "t2b-complex")
    {
        asymptoticLF<complexv> funct;
        funct.ReadF(fi1);
        funct.SaveF(fi2);
    }
    else if(key == "b2t-complex")
    {
        asymptoticLF<complexv> funct;
        funct.ReadF(fi1);
        funct.SaveFtxt(fi2);
    }
    else if(key == "t2b-complexdouble")
    {
        asymptoticLF_complexv funct;
        funct.ReadF(fi1);
        funct.SaveF(fi2);
    }
    else if(key == "b2t-complexdouble")
    {
        asymptoticLF_complexv funct;
        funct.ReadF(fi1);
        funct.SaveFtxt(fi2);
    }
    else
    {
        
        cout<<"Error: specified key is not recognized\n";
        
        cout<<"use: \n";
        cout<<"t2b-double\n";
        cout<<"b2t-double\n";
        cout<<"t2b-complex\n";
        cout<<"b2t-complex\n";
        cout<<"t2b-complexdouble\n";
        cout<<"b2t-complexdouble\n";
    }
    
     return 0;
}
