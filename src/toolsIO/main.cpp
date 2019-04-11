//
//  main.cpp
//  
//
//  Created by Darius Abramavicius on 5/22/12.
//  Copyright (c) 2012 VU. All rights reserved.
//

#include"../storage/storage.hpp"
#include"../dvector3d/dvector3d.hpp"
#include"../interpolationF/interpolationF.hpp"
#include"../asymptoticLF/asymptoticLF.hpp"


#include"toolsIO.hpp"





#include <iostream>
using std::cout;
#include <fstream>
using std::ifstream;
using std::ofstream;


// this program reads text from file
// and removes all commented lines.
// The comments are distinguished by # anywhere in the lines


int main(int argc, const char * argv[])
{
 //   // sample code
 //   toolsIO obj;
 //   
 //   cout<<obj.fromString<int>(string("12"))<<"\n";
    
    
    //  testing toWords
    //toolsIO* ttio = new toolsIO;
    //string istr = "  asdf tyty  -123456    +SDSDS56565   ";
    //int numw;
    //string* words = ttio->toWords(numw,istr);
    //cout<<"numw = "<<numw<<"\nwords:\n";
    //for(int ind=0; ind< numw; ind++)
    //    cout<<(words[ind]+"\n");
    //return 0;
    
    if(argc != 2)
    {
        cout<<"Error using toolsIO\n";
        return 0;
    }
    
    string inpf(argv[1]);
    
    ifstream* ifs = new ifstream(inpf.c_str());
        
    if(!ifs->is_open())
        cout<<"Error: toolsIO input file has not been openned\n";
    else
    {
        //cout<<"toolsIO: reading input file and removing comments\n";
            
            
        toolsIO* tio = new toolsIO;
        ofstream* ofs = new ofstream( (inpf+".nc").c_str());
        string str;
        
       
        while( !ifs->eof() )
        {
            str = "";
            tio->StreamSkipTrailers(ifs);
            getline(*ifs,str);
            tio->StreamSkipTrailers(str);
            *ofs<< str <<"\n";
        }

	// reading whole file into string
	std::stringstream strstr;
	tio->ReadWholeFile("toolsIO.hpp",strstr);
	cout<<strstr.str()<<"\n";

        delete tio;
        delete ofs;
    }
    delete ifs;
}
