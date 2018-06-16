
#include"calculator-main.hpp"

using std::cin;
using std::cout;

int main(int argc, const char * argv[])
{
    
    if(argc != 3 )
    {    cout<<"Error: specify the input file and the output file\n";
        return 0;
    }
    
    toolsIO tio;
    
    string ifilename(argv[1]);
    string ofilename(argv[2]);

    ifstream istr(ifilename.c_str());
    if(!istr.is_open())
    {
        cout<<"Error: input file cannot be read\nquitting.\n";
        return 0;
    }

    string calculator = "";
    tio.StreamSkipTrailers(istr);
    getline(istr,calculator);
    tio.StreamSkipTrailers(calculator);

    std::size_t found;
    found = calculator.find("3rd-2D-Reduced");
    if(found!=std::string::npos){
		// 2D calculators
		// perform calculations based on filename:
		calculator_3rd_reduced<propagatorExciton> calc(ifilename);
		calc.LaunchBasic();
		calc.Publish(ofilename);
		return 0;
	}
	
    found = calculator.find("3rd-2D-Secular");
    if(found!=std::string::npos){
		// 2D calculators
		// perform calculations based on filename:
		calculator_3rd_secular_cumulant calc(ifilename);
		calc.LaunchBasic();
		calc.Publish(ofilename);
		return 0;
	}

    found = calculator.find("3rd-2D-Wavefunctions");
    if(found!=std::string::npos){
		// 2D calculators
		// perform calculations based on filename:
		calculator_3rd_wf<propagatorExcitonWF,wavefunctionD1> calc(ifilename);
		calc.LaunchBasic();
		calc.Publish(ofilename);
		return 0;
	}
    found = calculator.find("1st-Absorption-Reduced");
    if(found!=std::string::npos){
		// absorption calculators

		// perform calculations based on filename:
		calculator_abs_reduced<propagatorExciton> calc(ifilename);
		calc.Launch();
		calc.Publish(ofilename);
		return 0;
	}

    found = calculator.find("1st-Absorption-Secular");
    if(found!=std::string::npos){
		// absorption calculators

		// perform calculations based on filename:
		calculator_abs_secular_cumulant calc(ifilename);
		calc.Launch();
		calc.Publish(ofilename);
		return 0;
	}

    found = calculator.find("1st-Absorption-Wavefunctions");
    if(found!=std::string::npos){
		// absorption calculators

		// perform calculations based on filename:
		calculator_abs_wf<propagatorExcitonWF,wavefunctionD1> calc(ifilename);
		calc.Launch();
		calc.Publish(ofilename);
		return 0;
	}
    

    
    
    
    return 0;
}
