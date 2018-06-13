#include<iostream>

#include"feinman2sideddiagram3.hpp"

#include"../toolsIO/toolsIO.hpp"

using std::cin;
using std::cout;

int main(int argc, const char * argv[])
{
    // this project calcuates the 
    // Feinman two-sided diagram at three time intervals
    
    // the diagram (from the header)
//       e |  | e
//        -------
//       k |  | l
//       d |  | f
//        -------
//       j |  | m
//       c |  | g
//        -------
//       i |  | n
//       b |  | h
//        -------
//       a |  | a

//                states[4 ] /*|*/ states[4 ] 
//                -------------------------
//                states[10] /*|*/ states[11] 
//                states[3 ] /*|*/ states[5 ] 
//                -------------------------
//                states[9 ] /*|*/ states[12] 
//                states[2 ] /*|*/ states[6 ] 
//                -------------------------
//                states[8 ] /*|*/ states[13] 
//                states[1 ] /*|*/ states[7 ] 
//                -------------------------
//                states[0 ] /*|*/ states[0 ] 

    // all is read from the input file
    if(argc != 3 )
    {    cout<<"Error: specify the input file and the output file\n";
        return 0;
    }

    toolsIO tio;
    string str = "";

    string ifilename(argv[1]);
    string ofilename(argv[2]);

    ifstream ifs(ifilename.c_str());

    
                
    cout<<"reading the t1 frequency (real part)\n";
    double freq1 = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                freq1  = tio.fromString<double>(str);

    cout<<"reading the t1 frequency (imag part)\n";
    double deph1 = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                deph1  = tio.fromString<double>(str);

    cout<<"reading the t2 frequency (real part)\n";
    double freq2 = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                freq1  = tio.fromString<double>(str);

    cout<<"reading the t2 frequency (imag part)\n";
    double deph2 = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                deph1  = tio.fromString<double>(str);

    cout<<"reading the t3 frequency (real part)\n";
    double freq3 = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                freq1  = tio.fromString<double>(str);

    cout<<"reading the t3 frequency (imag part)\n";
    double deph3 = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                deph1  = tio.fromString<double>(str);


    cout<<"reading four system transition dipoles of interactions\n";
    storage<dvector3d> dip(1); // two dimensional
                str = "";
                if(true){
                        storage<double> arr(2); // two dimensional
                        arr.Allocate(4,3);
                        if(tio.ReadRectangular(&arr, ifs))
                        {
                            cout<<"Error: Some problem with file reading\n";
                            return 1;
                        }
                        // making dipoles
                        dip.Allocate(4);
                        for(int ind=0; ind<4; ind ++)
                            dip.data1D[ind]=dvector3d(
                                    arr.data2D[ind][0],
                                    arr.data2D[ind][1],
                                    arr.data2D[ind][2]);
                        arr.Delete();
                }

    cout<<"reading four field  vectors of interactions\n";
    storage<dvector3d> fld(1); // two dimensional
                str = "";
                if(true){
                        storage<double> arr(2); // two dimensional
                        arr.Allocate(4,3);
                        if(tio.ReadRectangular(&arr, ifs))
                        {
                            cout<<"Error: Some problem with file reading\n";
                            return 1;
                        }
                        // making dipoles
                        fld.Allocate(4);
                        for(int ind=0; ind<4; ind ++)
                            fld.data1D[ind]=dvector3d(
                                    arr.data2D[ind][0],
                                    arr.data2D[ind][1],
                                    arr.data2D[ind][2]);
                        arr.Delete();
                }

                
                
                
    cout<<"reading the t1 time number of snapshots\n";
    int tim1N = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                tim1N  = tio.fromString<int>(str);

    cout<<"reading the t1 time final value\n";
    double tim1F = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                tim1F  = tio.fromString<double>(str);
                
                
                
    cout<<"reading the t2 time number of snapshots\n";
    int tim2N = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                tim2N  = tio.fromString<int>(str);

    cout<<"reading the t2 time final value\n";
    double tim2F = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                tim2F  = tio.fromString<double>(str);
                
                
                
    cout<<"reading the t3 time number of snapshots\n";
    int tim3N = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                tim3N  = tio.fromString<int>(str);

    cout<<"reading the t3 time final value\n";
    double tim3F = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                tim3F  = tio.fromString<double>(str);
                
                
                
                
                
    //////////////////
    //  optional fields
                
                
                               
    cout<<"checking if the diagram is coherent or transport: (0 coherent; 1 transport)\n";
    int transport = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                transport  = tio.fromString<int>(str);

    asymptoticLF<double>* transportGF = 0;
    if(transport)
    {
        // read the population transport Greens function
        cout<<"reading the population transport Greens function\n";
                // 1. allocation
                transportGF = new asymptoticLF<double>;

                // taking filename from the ifile
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                ifstream ifsl(str.c_str());
                transportGF->ReadF(&ifsl);
                ifsl.close();
    }
                
                
    // now the fluctuation part
    cout<<"checking if the cumulant fluctuations are included: (0 excluded; 1 included)\n";
    int cumulant = 0;
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                cumulant  = tio.fromString<int>(str);

    int numBO = 0;
    asymptoticLF_complexv* lineshapes  = 0;
    int numFL = 0;
    storage<double> ampFL(3); // THREE DIMENSIONAL
    int states[14];
    int technique;

    if(cumulant)
    {
        cout<<"taking the number of bath modes from the ifile\n";
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                numBO = tio.fromString<int>(str);
                //cout<<numBosc<<" test\n";

        cout<<"reading all lineshape functions from files\n";
                // 1. allocation
                delete[] lineshapes;
                lineshapes = new asymptoticLF_complexv[numBO];
                for(int ind = 0; ind<numBO; ind++)
                {
                        // taking filenames from the ifile
                        str = "";
                        tio.StreamSkipTrailers(&ifs);
                        getline(ifs,str);
                        tio.StreamSkipTrailers(str);
                        ifstream ifsl(str.c_str());
                        lineshapes[ind].ReadF(&ifsl);
                        ifsl.close();
                }
        
        cout<<"taking the technique definition index\n";
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                technique = tio.fromString<int>(str);
                //cout<<numBosc<<" test\n";


        cout<<"reading the interaction configuration (14 state indices on a single line)\n";
                str = "";
                if(true){
                        storage<double> arr(2); // two dimensional
                        arr.Allocate(1,14);
                        if(tio.ReadRectangular(&arr, ifs))
                        {
                            cout<<"Error: Some problem with file reading\n";
                            return 1;
                        }
                        for(int ind=0; ind<14; ind ++)
                            states[ind]=arr.data2D[0][ind];
                        arr.Delete();
                }
    

        cout<<"reading fluctuation amplitudes in accord with the configuration\n";
        cout<<"first: reading the number of elements in the fluctuation matrix\n";
                str = "";
                tio.StreamSkipTrailers(&ifs);
                getline(ifs,str);
                tio.StreamSkipTrailers(str);
                numFL  = tio.fromString<int>(str);

        cout<<"second: reading the  fluctuation matrices for each bath mode\n";
                str = "";
                ampFL.Allocate(numBO,numFL,numFL);
                for(int indo=0; indo<numBO; indo ++)
                {
                        storage<double> arr(2); // two dimensional
                        arr.Allocate(numFL,numFL);
                        if(tio.ReadRectangular(&arr, ifs))
                        {
                            cout<<"Error: Some problem with file reading\n";
                            return 1;
                        }
                        // saving  amplitudes
                        for(int indu=0; indu<numFL; indu ++)
                        for(int indl=0; indl<numFL; indl ++)
                            ampFL.data3D[indo][indu][indl]=arr.data2D[indu][indl];
                        arr.Delete();
                }

    }
                
    cout<<"Reading in finished\n";
                
    
    
    ifs.close();


    
    
    
    
    // making the diagram
    complexv ifreq[3];
    ifreq[0] = complexv(freq1, deph1);
    ifreq[1] = complexv(freq2, deph2);
    ifreq[2] = complexv(freq3, deph3);
    
    int averaging = 0;
    
    feinman2sideddiagram3 obj(ifreq, dip.data1D, fld.data1D,averaging);

    if(cumulant)
    {
	int first, second, third;
	ampFL.GetSize(first, second, third);
	int nstates = first + second;
        obj.assignLineshape(technique,lineshapes,ampFL.data3D,numBO,states, nstates);
    }

    if(transport)
        obj.assignPropagator(transportGF);

    // propagation from    0 to fin in each dimension 
    double tT[3];
    int tN[3];
    tT[0] = tim1F;
    tT[1] = tim2F;
    tT[2] = tim3F;
    tN[0] = tim1N;
    tN[1] = tim2N;
    tN[2] = tim3N;
    
    storage<complexv> res(3);
    res.Allocate(tim3N,tim2N,tim1N);
    
    
    cout<<"Performing propagation\n";
    
    obj.propagate( res, tT, tN);


    // writing the file
    ofstream ofs(ofilename.c_str());
    ofs<<"# feinman2sideddiagram3 project\n";
    ofs<<"# format:\n";
    ofs<<"# it3  it2  it1  res.real()  res.imag()\n";
    
    for(int it3 = 0; it3<tim3N; it3++)
    for(int it2 = 0; it2<tim2N; it2++)
    for(int it1 = 0; it1<tim1N; it1++)
    {
        ofs<<it3<<"\t"<<it2<<"\t"<<it1<<"\t";
        ofs<<res.data3D[it3][it2][it1].real()<<"\t";
        ofs<<res.data3D[it3][it2][it1].imag()<<"\n";
    }
    
    ofs.close();     

    
    
    // deleting additional allocations
    if(transport)
        delete transportGF;
    
    if(cumulant)
        delete[] lineshapes;


    
    return 0;
}
