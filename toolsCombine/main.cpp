//
//  main.cpp
//  
//
//  Created by Darius Abramavicius on 5/22/12.
//  Copyright (c) 2012 VU. All rights reserved.
//

#include"toolsCombine.hpp"
#include"../toolsIO/toolsIO.hpp"
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
    //combining two files of absorption
    
    toolsIO tio;
    ifstream ifs;
    ofstream ofs;

    if(argc != 5 )
    {    cout<<"Error: specify the number of values, input file 1, input file 2 and the output file\n";
        return 0;
    }
 
    
    int num = tio.fromString<int>(string(argv[1]));
    string ifilename1(argv[2]);
    string ifilename2(argv[3]);
    string ofilename(argv[4]);
    
                storage<double> dat1(2);
                dat1.Allocate(num,2);
    
                ifs.open(ifilename1.c_str());
                if( tio.ReadRectangular(&dat1, ifs))
                {
                    cout<<"Error: Some problem with file reading\n";
                    return 0;// 1;
                }
                ifs.close();
                
                //for(int ind = 0; ind<num; ind++)
                //   cout<<dat1.data2D[ind][1]<<"\n";


                storage<double> dat2(2);
                dat2.Allocate(num,2);
    
                ifs.open(ifilename2.c_str());
                if( tio.ReadRectangular(&dat2, ifs))
                {
                    cout<<"Error: Some problem with file reading\n";
                    return 0;// 1;
                }
                ifs.close();
    
                
                // adding the data
                for(int ind = 0; ind<num; ind++)
                    dat2.data2D[ind][1] += dat1.data2D[ind][1];
                            
                // printing the result
                cout<<"writing the result\n";
                ofs.open(ofilename.c_str());
                tio.WriteRectangular(&dat2, ofs);
                ofs.close();

}
