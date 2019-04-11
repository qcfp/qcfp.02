#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include"../storage/storage.hpp"

using std::ifstream;
using std::istream;
using std::ofstream;
using std::string;
using std::istringstream;
using std::ostringstream;




class toolsIO
{
public:
    
	// prepares stream for the meaningfull data input
    // by removing front delimiters
    // and end comments '#'
	int StreamSkipTrailers(istream& str);
	int StreamSkipTrailers(ifstream* str);
    int StreamSkipTrailers(string& sstr);
    int StreamSkipTrailers(istringstream& sstr);
    int StreamSkipTrailers(istringstream* sstr);

    template <class T>
    static T fromString(const string& s)
    {
        istringstream stream(s);
        T t;
        stream >> t;
        return t;
    }
    template<class T>
    static string toString(const T& t)
    {
        ostringstream stream;
        stream << t;
        return stream.str();
    }


	int ReadWholeFile(std::string name,std::stringstream& strstr);
    
    // next function removes commented lines
	int ReadWholeFile(std::string name,std::stringstream& strstr,int flag);

    
    // next takes in the string of words
    // and splits it into separate words  given as a dynamic array
    // based on "space" and "tab" delimiters
    string*  toWords(int& num, string& istr);
    // on output function returns the array of meaningful  words
    // (do not forget to delete it)
    // "num" is an output number of words
    // "istr" is an input string
    
    template<class T>
    int LookUpAndReadSquareMatrix(string keyword,
                                  string label, 
                                  int numl, int numr, 
                                  storage<T>& arr,
                                string ifile)
    {
        if(arr.CheckDimension() != 2){
            cout<<"Error: reading into wrong dimensional matrix: skipping\n";
            return 0;
        }
        string& sstr = ifile;
        string str = "";
        if(arr.IsAlloc())
            arr.Delete();

        // looking for entry
        std::size_t found = sstr.find(keyword);
        if(found != std::string::npos)
        {
            sstr = sstr.substr(found);
            std::istringstream ifs(sstr);
            getline(ifs,keyword);
            cout<<label;

            arr.Allocate(numl,numr);
            if(ReadRectangular(arr, ifs))
            {
                cout<<"Error: Some problem with file reading\n";
                return 0;// 1;
            }
            return 1;
        }
        return 0;
    }


    template<class T>
    int LookUpAndReadTriangularMatrix(string keyword,
                                  string label, 
                                  int numl,
                                  storage<T>& arr,
                                string ifile)
    {
        if(arr.CheckDimension() != 2){
            cout<<"Error: reading into wrong dimensional matrix: skipping\n";
            return 0;
        }
        string& sstr = ifile;
        string str = "";

        if(arr.IsAlloc())
            arr.Delete();

        int& numr=numl;

        // looking for entry
        std::size_t found = sstr.find(keyword);
        if(found != std::string::npos)
        {
            sstr = sstr.substr(found);
            std::istringstream ifs(sstr);
            getline(ifs,keyword);
            cout<<label;

            arr.Allocate(numl,numr);
            if(ReadTriangular(arr, ifs))
            {
                  cout<<"Error: Some problem with file reading\n";
                  return 0;// 1;
            }
            return 1;
        }
        return 0;
    }

   
    
    // a Function to read a triangular dataset from a file
    // with predefined size
    int ReadTriangular(storage<double>* result, string ifi);

    // a Function to read a square dataset from a file
    // with predefined size
    int ReadRectangular(storage<double>* result, string ifi);
    
    // a Function to write a triangular dataset to a file
    // with predefined size
    int WriteTriangular(storage<double>* result, string ifi);
    
    // a Function to write a square dataset to a file
    // with predefined size
    int WriteRectangular(storage<double>* result, string ifi);
    
//----------------
    
    // a Function to read a triangular dataset from a file
    // with predefined size
    int ReadTriangular(storage<double>* result, ifstream& ifi);
    int ReadTriangular(storage<double>& result, istream& ifi);

    // a Function to read a square dataset from a file
    // with predefined size
    int ReadRectangular(storage<double>* result, ifstream& ifi);
    int ReadRectangular(storage<double>& result, istream& ifi);
    
    // a Function to write a triangular dataset to a file
    // with predefined size
    int WriteTriangular(storage<double>* result, ofstream& ifi);
    int WriteTriangular(storage<double>& result, ofstream& ifi);
    
    // a Function to write a square dataset to a file
    // with predefined size
    int WriteRectangular(storage<double>* result, ofstream& ifi);
    int WriteRectangular(storage<double>& result, ofstream& ifi);

//------------------
    // a Function to read a triangular dataset from a file
    // with predefined size
    int ReadTriangular(storage<int>& result, istream& ifi);

    // a Function to read a square dataset from a file
    // with predefined size
    int ReadRectangular(storage<int>& result, istream& ifi);
    
    // a Function to write a triangular dataset to a file
    // with predefined size
    int WriteTriangular(storage<int>& result, ofstream& ifi);
    
    // a Function to write a square dataset to a file
    // with predefined size
    int WriteRectangular(storage<int>& result, ofstream& ifi);
    
private:
    int llen(char*,int);
    
};
