#include"../storage/storage.hpp"

#include"toolsIO.hpp"

#include <string>
#include <fstream>
#include <ostream>
#include <iostream>

#include<cctype>
#include <ctype.h>

using std::ifstream;
using std::ofstream;
using std::string;
using std::istringstream;
using std::cout;

int toolsIO::ReadWholeFile(std::string name,std::stringstream& strstr,int flag)
{
    std::stringstream temporary;
    std::string str;
    ReadWholeFile(name,temporary);
        if(flag) 
        {
            // and now remove commented lines
            while(!temporary.eof())
            {
                    StreamSkipTrailers(temporary);
                    getline(temporary,str);
                    //StreamSkipTrailers(str);
                    strstr<<str<<"\n";                    
            }
        }
}

int toolsIO::ReadWholeFile(std::string name,std::stringstream& strstr)
{
	std::ifstream istr(name.c_str(),std::ios::binary);
	//std::streambuf* puf;

	if(istr.is_open())
	{
		strstr<<istr.rdbuf();
		return 0;
	}
	std::cout<<"Error: file reading in toolsIO::ReadWholeFile failed\n";
	return 1;
}
int toolsIO::StreamSkipTrailers(ifstream* str)
{
	return StreamSkipTrailers((istream&) *str);
}
int toolsIO::StreamSkipTrailers(istream& str)
{
    char stemp[65536]={0};
    std::string tstr;
    int err = 0;
    while(!str.eof())
    {// skipping comments
        stemp[1] = 0;
        stemp[0] = str.peek();// checking first character
        if(stemp[0]=='#')
        {
            getline(str,tstr);
            //do // removing any occurrence of loooong paragraph
            //{
            //    str->getline(stemp,65535);
            //}while(llen(stemp,65536) == 65535 );
        }
        else if(stemp[0]=='\n')// skipping empty lines
            //str->getline(stemp,65535);
            str.get();
        else if(stemp[0]=='\r')// skipping tariling spaces
            str.get();
        else if(stemp[0]==' ')// skipping tariling spaces
            str.get();
        else if(stemp[0]=='\t')// skipping tariling spaces
            str.get();
        else if(std::isprint(stemp[0]) )
            break;
        else
        {
            str.get();
            // some unknown character
            err ++;
            if(err>1000)
            {
                std::cout<<"Error in file leading: too many unrecognized characters\n";
                break;
            }
        }
    }
    return 0;
}
\
int toolsIO::StreamSkipTrailers(istringstream& str)
{
	return StreamSkipTrailers(&str);
}
int toolsIO::StreamSkipTrailers(istringstream* str)
{
    string tstr = str->str();
    int numc = tstr.length();
    int numext = 0;
    //cout<<"StreamSkipTrailers: "<<numc<<"\n";
    char stemp[2]={0};
    while(true)
    {// skipping comments
        stemp[1] = 0;
        stemp[0] = str->peek();// checking first character
        if(stemp[0]=='#')
        {
           // removing all the rest
            str->str(string(""));
            break;
        }
        //else if(stemp[0]=='\n')// skipping empty lines
        //    str->getline(stemp,65535);
        else if(stemp[0]==' ')// skipping tariling spaces
        {
            str->get();
            numext ++;
        }
        else if(stemp[0]=='\t')// skipping tariling spaces
        {
            str->get();
            numext ++;
        }
        else if(stemp[0]=='\r')// skipping tariling spaces
        {
            str->get();
            numext ++;
        }
        else if(!std::isprint(stemp[0]) )
        {
            str->get();
            numext ++;
        }
        else break;

    
        if(numext == numc)
        {
            // extracted all characters
            str->str(string(""));
            break;
        }
    }
    return 0;
}


int toolsIO::StreamSkipTrailers(string& sstr)
{
    //istringstream str(sstr);
    // removing frontal spaces
    //StreamSkipTrailers(&str);


    
    std::size_t found; 
    

    // removing comments
    found = sstr.find_first_of("#");
    if(found != std::string::npos)
        sstr = sstr.substr(0,found);

    //  scanning of any spaces or special characters at the beginning
    while(true)
    {
        found = sstr.find_first_of(" ");
        if(found != std::string::npos)
            sstr = sstr.substr(found+1);
        else break;
    }

    while(true)
    {
        found = sstr.find_first_of("\t");
        if(found != std::string::npos)
            sstr = sstr.substr(found+1);
        else break;
    }    
    
    
    // finally scanning of any special characters at the end
    found = sstr.find_first_of("\r");
    if(found != std::string::npos)
        sstr = sstr.substr(0,found);
    
    
    //std::cout <<sstr << '\n';
        
    return 0;
}

int toolsIO::llen(char* ia, int in)
{
    int len=0;
    while ( ia[len]!= '\0' ) {
        len++;
        if(len==in) break;
    }
    return len;
}

string* toolsIO::toWords(int& num, string& istr)
{
    ////// this code is to make words from string istr
    ////// "words" must be unallocated
    
    // creating output buffer
    int mult = 1;
    string* words = new string[100];// first guessing no more than 100 words
    num = 0;// number of currently found words
    
    //studying the input string
    string& singleLine = istr;
    int length = singleLine.length(); 
    //cout<<singleLine<<"\n";
    
    // removing front delimiters
    StreamSkipTrailers(singleLine);
    length = singleLine.length(); 
    //cout<<singleLine<<"\n";
    
    // nonzero string
    do
    {
        if(length == 0) return words;
        
        // searching for a delimiter
        size_t f1 = singleLine.find(" ");
        size_t f2 = singleLine.find("\t");
        size_t found = f2; if(f1<f2) found = f1;
        
        if(found == string::npos)
        {
            //no more delimiters 
            words[num]=singleLine.substr(0,length);
            num ++;
            return words;
        }
        else
        {
            // getting a word and moving forward
            words[num]=singleLine.substr(0,found);
            num ++;
            singleLine=singleLine.substr(found);
            StreamSkipTrailers(singleLine);
            length = singleLine.length();
            //cout<<length<<"\n";
        }
        
        
        // analyzing the number of words found
        if(num== mult*100)
        {
            // have to increase the length of words
            mult ++;
            
            //if(mult == 100) // too many words
            {
                cout<<"Error: toolsIO found too many words\n";
                return 0;
            }
            
            string* twords = new string[mult*100];
            for(int ti=0; ti<num; ti++)
            {
                twords[ti]=words[ti];
            }
            delete[] words;
            words=twords;
        }
        
        
    }while(true);
    
    

}


// a Function to read a triangular dataset from a file
// with predefined size
int toolsIO::ReadTriangular(storage<double>* result, string ifi)
// result should be allocated
{
   	if (ifi == "" )
    {
        cout<<"Warning: request to write to \"no-name\" file: skipping\n";
        return 0;
    }
    if(result == 0)
    {
        cout<<"Error: request to place data in unallocated space.\n";
        return 1;
    }
    int dim1,dim2;
    result->GetSize(dim1,dim2);
    if(dim1==0 || dim2==0)
    {
        cout<<"Warning: skipping reading of zero-size data.\n";
        return 0;
    }
    ifstream* ifi_str = new ifstream(ifi.c_str());
    ReadTriangular(result, *ifi_str);
    
    ifi_str->close();
    delete ifi_str;
    
    return 0;
    
}
// a Function to read a square dataset from a file
// with predefined size
int toolsIO::ReadRectangular(storage<double>* result, string ifi)
// result should be allocated
{
   	if (ifi == "" )
    {
        cout<<"Warning: request to write to \"no-name\" file: skipping\n";
        return 0;
    }
    if(result == 0)
    {
        cout<<"Error: request to place data in unallocated space.\n";
        return 1;
    }
    
    
    if(result->CheckDimension() == 1)
    {
        int dim1,dim2;
        result->GetSize(dim1);
        if(dim1==0)
        {
            cout<<"Warning: skipping reading of zero-size data.\n";
            return 0;
        }
        
    }
    else if(result->CheckDimension() == 2)
    {
        int dim1,dim2;
        result->GetSize(dim1,dim2);
        if(dim1==0 || dim2==0)
        {
            cout<<"Warning: skipping reading of zero-size data.\n";
            return 0;
        }
        
    }
    
    
    
 	ifstream* ifi_str = new ifstream(ifi.c_str());
        ReadRectangular(result, *ifi_str);
    
        ifi_str->close();
        delete ifi_str;
        
        return 0;
    
}

int toolsIO::WriteRectangular(storage<double>* result, string ifi)
// result should be allocated
{
   	if (ifi == "" )
    {
        cout<<"Warning: request to write to \"no-name\" file: skipping\n";
        return 1;
    }
    if(result == 0)
    {
        cout<<"Error: request to write unallocated data.\n";
        return 1;
    }
    
    
    if(result->CheckDimension() == 1)
    {
        int dim1,dim2;
        result->GetSize(dim1);
        if(dim1==0)
        {
            cout<<"Warning: skipping reading of zero-size data.\n";
            return 0;
        }
        
    }
    else if(result->CheckDimension() == 2)
    {
        int dim1,dim2;
        result->GetSize(dim1,dim2);
        if(dim1==0 || dim2==0)
        {
            cout<<"Warning: skipping reading of zero-size data.\n";
            return 0;
        }
        
    }

    ofstream* ifi_str = new ofstream(ifi.c_str());
    WriteRectangular(result, *ifi_str);
    
    ifi_str->close();
    delete ifi_str;
    
    return 0;
    
}

int toolsIO::WriteTriangular(storage<double>* result, string ifi)
// result should be allocated
{
   	if (ifi == "" )
    {
        cout<<"Warning: request to write to \"no-name\" file: skipping\n";
        return 1;
    }
    if(result == 0)
    {
        cout<<"Error: request to write unallocated data.\n";
        return 1;
    }
    int dim1,dim2;
    result->GetSize(dim1,dim2);
    if(dim1==0 || dim2==0)
    {
        cout<<"Warning: skipping writing of zero-size data.\n";
        return 0;
    }
    ofstream* ifi_str = new ofstream(ifi.c_str());
    WriteTriangular(result, *ifi_str);
    
    ifi_str->close();
    delete ifi_str;
    
    return 0;
    
}



// a Function to read a triangular dataset from a file
// with predefined size
int toolsIO::ReadTriangular(storage<double>* result, ifstream& ifi_str)
// result should be allocated
{
    if(!ifi_str.is_open())
	{
		cout<<"Error: triangular matrix file was not openned\n";
        return 1;
	}
	return ReadTriangular(*result, (istream&) ifi_str);
}    


// a Function to read a triangular dataset from a file
// with predefined size
int toolsIO::ReadTriangular(storage<double>& result, istream& ifi_str)
// result should be allocated
{
    
    int numv=1;
    if(!result.IsAlloc())
    {
        cout<<"Warning: request to read into unallocated matrix: skipping\n";
        return 1;
    }
    if(result.CheckDimension() != 2)
    {
        cout<<"Error: request to read into wrong dimensionality matrix: skipping\n";
        return 1;
    }
    
    result.GetSize(numv,numv);
    
    if(numv == 0)
    {
        cout<<"Warning: request to read zero size matrix\n";
        return 0;
    }
    
    cout<<"Reading matrix "<<numv<<" x "<<numv<<"\n";
    for(int inda = 0; inda<  numv; inda++)
		for(int indb = 0; indb<= inda; indb++)
        {
            StreamSkipTrailers(ifi_str);
            ifi_str >> result.data2D[inda][indb];
                cout<<result.data2D[inda][indb];
                if(indb == inda) cout<<"\n"; else cout<<"\t";
            //cout<<"read "<<result->data2D[inda][indb]<<"\n";
        }
    return 0;
    
}
// a Function to read a square dataset from a file
// with predefined size
int toolsIO::ReadRectangular(storage<double>* result, ifstream& ifi_str)
// result should be allocated
{
    
    
    if(!ifi_str.is_open())
    {
        cout<<"Error: triangular matrix file was not openned\n";
        return 1;
    }
    ReadRectangular(*result, (istream&) ifi_str);
}
    
int toolsIO::ReadRectangular(storage<double>& result, istream& ifi_str)
// result should be allocated
{
    
    
    int numvl=1;
    int numvr=1;
    int type;
    if(!result.IsAlloc())
    {
        cout<<"Warning: request to read into unallocated matrix: skipping\n";
        return 1;
    }
    

    
    if(result.CheckDimension() == 1)
    {
        type = 1;
        result.GetSize(numvl);
    }
    else if(result.CheckDimension() == 2)
    {
        type = 2;
        result.GetSize(numvl,numvr);
    }
    else
    {
        cout<<"Warning: request to write wrong dimension matrix: skipping\n";
        return 1;
    }
    
    
    
    if(numvl == 0 || numvr == 0 )
    {
        cout<<"Warning: request to read zero size matrix\n";
        return 0;
    }
    
    if(type == 2)
    {
        cout<<"Reading matrix "<<numvl<<" x "<<numvr<<"\n";
        for(int inda = 0; inda<  numvl; inda++)
            for(int indb = 0; indb< numvr; indb++)
            {
                StreamSkipTrailers(ifi_str);
                ifi_str >> result.data2D[inda][indb];
                cout<<result.data2D[inda][indb];
                if(indb == numvr-1) cout<<"\n"; else cout<<"\t";
            }
        return 0;
    }
    if(type == 1)
    {
        cout<<"Reading vector "<<numvl<<"\n";
        for(int inda = 0; inda<  numvl; inda++)
        {
            StreamSkipTrailers(ifi_str);
            ifi_str >> result.data1D[inda];
            cout<<result.data1D[inda];
            cout<<"\n";
        }
        return 0;
    }
    
    return 0;
}

int toolsIO::WriteRectangular(storage<double>* result, ofstream& ifi_str)
// result should be allocated
{
    int numvl=1;
    int numvr=1;
    int type;

    
    
    if(!ifi_str.is_open())
    {
        cout<<"Error: matrix file was not openned for writing\n";
        return 1;
    }
    
    
    if(result == 0)
    {
        cout<<"Error: request to write unallocated data.\n";
        return 1;
    }
    
    if(!result->IsAlloc())
    {
        cout<<"Warning: request to write unallocated matrix: skipping\n";
        return 1;
    }

    if(result->CheckDimension() == 1)
    {
        type = 1;
        result->GetSize(numvl);
    }
    else if(result->CheckDimension() == 2)
    {
        type = 2;
        result->GetSize(numvl,numvr);
    }
    else
    {
        cout<<"Warning: request to write wrong dimension matrix: skipping\n";
        return 1;
    }
    
    
    if(numvl == 0 || numvr == 0 )
    {
        cout<<"Warning: request to write zero size matrix: skipping\n";
        return 1;
    }
    
    
    if(type == 2)
    {
        for(int inda = 0; inda<  numvl; inda++)
        {
            for(int indb = 0; indb< numvr-1; indb++)
            {
                ifi_str << result->data2D[inda][indb]<<"\t";
            }
            ifi_str << result->data2D[inda][numvr-1]<<"\n";
        }
    }

    if(type == 1)
    {
        for(int inda = 0; inda<  numvl; inda++)
        {
            ifi_str << result->data1D[inda]<<"\n";
        }
    }
        
    return 0;
    
}

int toolsIO::WriteTriangular(storage<double>* result, ofstream& ifi_str)
// result should be allocated
{
    int numvl=1;
    int numvr=1;
    
    if(result == 0)
    {
        cout<<"Error: request to write unallocated data.\n";
        return 1;
    }
    if(!result->IsAlloc())
    {
        cout<<"Warning: request to write unallocated matrix: skipping\n";
        return 1;
    }

    if(result->CheckDimension() == 2)
    {
        result->GetSize(numvl,numvr);
    }
    else
    {
        cout<<"Error: request to write wrong dimension matrix: skipping\n";
        return 1;
    }
    
    
    if(numvl == 0 || numvr != numvl )
    {
        cout<<"Error: request to write zero or non-square size matrix: skipping\n";
        return 1;
    }
    
    
        
    if(!ifi_str.is_open())
	{
		cout<<"Error: matrix file was not openned for writing\n";
        return 1;
	}
    
    for(int inda = 0; inda<  numvl; inda++)
    {
        for(int indb = 0; indb< inda; indb++)
        {
            ifi_str << result->data2D[inda][indb]<<"\t";
        }
        ifi_str << result->data2D[inda][inda]<<"\n";
    }
    
    return 0;
    
}




//----------------------------
// a Function to read a triangular dataset from a file
// with predefined size
int toolsIO::ReadTriangular(storage<int>& result, istream& ifi_str)
// result should be allocated
{
    
    int numv=1;
    if(!result.IsAlloc())
    {
        cout<<"Warning: request to read into unallocated matrix: skipping\n";
        return 1;
    }
    if(result.CheckDimension() != 2)
    {
        cout<<"Error: request to read into wrong dimensionality matrix: skipping\n";
        return 1;
    }
    
    result.GetSize(numv,numv);
    
    if(numv == 0)
    {
        cout<<"Warning: request to read zero size matrix: skipping\n";
        return 0;
    }
    
    cout<<"Reading matrix "<<numv<<" x "<<numv<<"\n";
    for(int inda = 0; inda<  numv; inda++)
        for(int indb = 0; indb<= inda; indb++)
        {
            StreamSkipTrailers(ifi_str);
            ifi_str >> result.data2D[inda][indb];
            cout<<result.data2D[inda][indb];
            if(indb == inda) cout<<"\n"; else cout<<"\t";
            //cout<<"read "<<result->data2D[inda][indb]<<"\n";
        }
        return 0;
    
}
    
int toolsIO::ReadRectangular(storage<int>& result, istream& ifi_str)
// result should be allocated
{
    
    
    int numvl=1;
    int numvr=1;
    int type;

    if(!result.IsAlloc())
    {
        cout<<"Warning: request to read into unallocated matrix: skipping\n";
        return 1;
    }
    
    
    if(result.CheckDimension() == 1)
    {
        type = 1;
        result.GetSize(numvl);
    }
    else if(result.CheckDimension() == 2)
    {
        type = 2;
        result.GetSize(numvl,numvr);
    }
    else
    {
        cout<<"Warning: request to write wrong dimension matrix: skipping\n";
        return 1;
    }
    
    
    
    if(numvl == 0 || numvr == 0 )
    {
        cout<<"Warning: request to read zero size matrix\n";
        return 0;
    }
    
    if(type == 2)
    {
        cout<<"Reading matrix "<<numvl<<" x "<<numvr<<"\n";
        for(int inda = 0; inda<  numvl; inda++)
            for(int indb = 0; indb< numvr; indb++)
            {
                StreamSkipTrailers(ifi_str);
                ifi_str >> result.data2D[inda][indb];
                cout<<result.data2D[inda][indb];
                if(indb == numvr-1) cout<<"\n"; else cout<<"\t";
            }
        return 0;
    }
    if(type == 1)
    {
        cout<<"Reading vector "<<numvl<<"\n";
        for(int inda = 0; inda<  numvl; inda++)
        {
            StreamSkipTrailers(ifi_str);
            ifi_str >> result.data1D[inda];
            cout<<result.data1D[inda];
            cout<<"\n";
        }
        return 0;
    }
    
    return 0;
}

int toolsIO::WriteRectangular(storage<int>& result, ofstream& ifi_str)
// result should be allocated
{
    int numvl=1;
    int numvr=1;
    int type;

    if(!result.IsAlloc())
    {
        cout<<"Warning: request to write unallocated matrix: skipping\n";
        return 1;
    }
    
    if(!ifi_str.is_open())
    {
        cout<<"Error: matrix file was not openned for writing\n";
        return 1;
    }
    
    
    if(result.CheckDimension() == 1)
    {
        type = 1;
        result.GetSize(numvl);
    }
    else if(result.CheckDimension() == 2)
    {
        type = 2;
        result.GetSize(numvl,numvr);
    }
    else
    {
        cout<<"Warning: request to write wrong dimension matrix: skipping\n";
        return 1;
    }
    
    
    if(numvl == 0 || numvr == 0 )
    {
        cout<<"Warning: request to write zero size matrix: skipping\n";
        return 1;
    }
    
    
    if(type == 2)
    {
        for(int inda = 0; inda<  numvl; inda++)
        {
            for(int indb = 0; indb< numvr-1; indb++)
            {
                ifi_str << result.data2D[inda][indb]<<"\t";
            }
            ifi_str << result.data2D[inda][numvr-1]<<"\n";
        }
    }

    if(type == 1)
    {
        for(int inda = 0; inda<  numvl; inda++)
        {
            ifi_str << result.data1D[inda]<<"\n";
        }
    }
        
    return 0;
    
}

int toolsIO::WriteTriangular(storage<int>& result, ofstream& ifi_str)
// result should be allocated
{
    int numvl=1;
    int numvr=1;
    
    if(!result.IsAlloc())
    {
        cout<<"Warning: request to write unallocated matrix: skipping\n";
        return 1;
    }
    
    if(result.CheckDimension() == 2)
    {
        result.GetSize(numvl,numvr);
    }
    else
    {
        cout<<"Error: request to write wrong dimension matrix: skipping\n";
        return 1;
    }
    
    
    if(numvl == 0 || numvr != numvl )
    {
        cout<<"Error: request to write zero or non-square size matrix: skipping\n";
        return 1;
    }
    
    
        
    if(!ifi_str.is_open())
	{
		cout<<"Error: matrix file was not openned for writing\n";
        return 1;
	}
    
    for(int inda = 0; inda<  numvl; inda++)
    {
        for(int indb = 0; indb< inda; indb++)
        {
            ifi_str << result.data2D[inda][indb]<<"\t";
        }
        ifi_str << result.data2D[inda][inda]<<"\n";
    }
    
    return 0;
    
}
