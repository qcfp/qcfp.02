#include<iostream>
#include<fstream>

#include"constructor-disorder.hpp"
#include"../toolsIO/toolsIO.hpp"

using std::cin;
using std::cout;
using std::ifstream;
using std::ofstream;

int main(int argc, const char * argv[])
{
    toolsIO tio;

    if(argc != 7 )
    {    cout<<"Error: specify the seed, the key-prefix, the key-suffix, the number of elements, the input file and the output file\n";
        return 0;
    }
 
    
    int iseed = tio.fromString<int>(string(argv[1]));
    int ikeyp = tio.fromString<int>(string(argv[2]));
    string ikeys(argv[3]); 
    int numEl = tio.fromString<int>(string(argv[4]));
    string ifilename(argv[5]);
    string ofilename(argv[6]);
    
    
    if(true)// diagonal disorder to the triangular matrix
    {
                storage<double> ham(2);
    
                cout<<"reading the triangular Matrix\n";
                ham.Allocate(numEl,numEl);

                
                ifstream ifs(ifilename.c_str());
                if( tio.ReadTriangular(&ham, ifs))
                {
                    cout<<"Error: Some problem with file reading\n";
                    return 0;// 1;
                }
                ifs.close();
                
                
                cout<<"Adding the diagonal disorder\n";
                constructor_disorder cds(iseed);
                cds.DiagIntra(ham, tio.fromString<double>(ikeys));

                cout<<"writing the triangular Matrix\n";
                ofstream ofs(ofilename.c_str());
                tio.WriteTriangular(&ham, ofs);
                ofs.close();
    }
 
    return 0;
}
