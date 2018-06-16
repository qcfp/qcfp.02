#include"toolsFFT.hpp"
#include"../constants/constants.hpp"


toolsFFT::toolsFFT()
{
	nump = 0;
	nump2 = 0;
    dimension = 0;

	FFTdatI = 0;
	FFTdatO = 0;

	fplanP = 0;
	fplanN = 0;

	locked = 0;

}
toolsFFT::~toolsFFT()
{
	cleanFFTSpace();
}


void toolsFFT::prepareFFTSpace(int np)
{
	nump = np;
    dimension = 1;
    
    // onedimenional transforms
    FFTdatI = (fftw_complex*)fftw_malloc(nump*sizeof(fftw_complex));
    FFTdatO = (fftw_complex*)fftw_malloc(nump*sizeof(fftw_complex));
    fplanP=fftw_plan_dft_1d(nump,FFTdatI,FFTdatO,FFTW_FORWARD,FFTW_ESTIMATE);
    fplanN=fftw_plan_dft_1d(nump,FFTdatI,FFTdatO,FFTW_BACKWARD,FFTW_ESTIMATE);
    locked = 1;
}
void toolsFFT::prepareFFTSpace(int in1, int in2)
{
	nump = in1;
    nump2 = in2;
    dimension = 2;
    
    // preparing for 2D Fourier transforms
    FFTdatI = (fftw_complex*)fftw_malloc(in1*in2*sizeof(fftw_complex));
    FFTdatO = (fftw_complex*)fftw_malloc(in1*in2*sizeof(fftw_complex));
    fplanP=fftw_plan_dft_2d(in1,in2,FFTdatI,FFTdatO,FFTW_FORWARD,FFTW_ESTIMATE);
    fplanN=fftw_plan_dft_2d(in1,in2,FFTdatI,FFTdatO,FFTW_BACKWARD,FFTW_ESTIMATE);
    locked = 1;
}

void toolsFFT::cleanFFTSpace()
{

    nump = 0;
    nump2 = 0;
    locked = 0;
    dimension = 0;
    fftw_free(FFTdatI); 
    fftw_free(FFTdatO);
    fftw_destroy_plan(fplanP);
    fftw_destroy_plan(fplanN);
}

void toolsFFT::executeGen(int& selection)
{
	if(selection == 1)
		fftw_execute(fplanP);
	else if(selection == 0)
		fftw_execute(fplanN);
}

void toolsFFT::executeP(complexv* result,complexv* source, int np)
//function for 1D requested
{
	int selection = 1; // forward transform requested

	if(locked  == 0)
		// preparations
		prepareFFTSpace(np);
        // now it is prepared
    
	if( nump == np) // OK proceeding
	{
		// copying data
		for(int in = 0; in <nump; in ++)
		{
			FFTdatI[in][0]=source[in].real();
			FFTdatI[in][1]=source[in].imag();
                        
                        //cout<<FFTdatI[in][0]<<" fft I\n";
                        //cout<<FFTdatI[in][1]<<" fft I\n";
		}
		fftw_execute(fplanP);

		// saving the calculated result
		for(int indt = 0; indt < nump; indt ++)
                {
			result[indt] = complexv( FFTdatO[indt][0], FFTdatO[indt][1] );
                        //cout<<FFTdatO[indt][0]<<" fft O\n";
                        //cout<<FFTdatO[indt][1]<<" fft O\n";
                }
	}
	else
	{
		cleanFFTSpace();
		executeP(result,source,np);
	}
}
void toolsFFT::executeP(complexv** result,complexv** source, int in1, int in2)
//function for 2D requested
{
	int selection = 1; // forward transform requested
    
	if(locked  == 0)
		// preparations
		prepareFFTSpace(in1,in2);
        // now it is prepared
    
	if( nump == in1 && nump2 == in2) // OK proceeding
	{
		// copying data
		for(int in1 = 0; in1 <nump; in1 ++)
            for(int in2 = 0; in2 <nump2; in2 ++)
		{
			FFTdatI[in1*nump2+in2][0]=source[in1][in2].real();
			FFTdatI[in1*nump2+in2][1]=source[in1][in2].imag();
		}
		fftw_execute(fplanP);
        
		// saving the calculated result
		for(int in1 = 0; in1 <nump; in1 ++)
            for(int in2 = 0; in2 <nump2; in2 ++)
			result[in1][in2] = complexv( FFTdatO[in1*nump2+in2][0], FFTdatO[in1*nump2+in2][1] );
	}
	else
	{
		cleanFFTSpace();
		executeP(result,source,in1,in2);
	}
}


void toolsFFT::executeN(complexv* result,complexv* source, int np)
//function for 1D requested
{
	int selection = 0; // back transform requested
    
	if(locked  == 0)
		// preparations
		prepareFFTSpace(np);
    // now it is prepared
    
	if( nump == np) // OK proceeding
	{
		// copying data
		for(int in = 0; in <nump; in ++)
		{
			FFTdatI[in][0]=source[in].real();
			FFTdatI[in][1]=source[in].imag();
		}
		fftw_execute(fplanN);
        
		// saving the calculated result
		for(int indt = 0; indt < nump; indt ++)
			result[indt] = complexv( FFTdatO[indt][0], FFTdatO[indt][1] );
	}
	else
	{
		cleanFFTSpace();
		executeN(result,source,np);
	}
}

void toolsFFT::executeN(complexv** result,complexv** source, int in1, int in2)
//function for 2D requested
{
	int selection = 0; // forward transform requested
    
	if(locked  == 0)
		// preparations
		prepareFFTSpace(in1,in2);
    // now it is prepared
    
	if( nump == in1 && nump2 == in2) // OK proceeding
	{
		// copying data
		for(int in1 = 0; in1 <nump; in1 ++)
            for(int in2 = 0; in2 <nump2; in2 ++)
            {
                FFTdatI[in1*nump2+in2][0]=source[in1][in2].real();
                FFTdatI[in1*nump2+in2][1]=source[in1][in2].imag();
            }
		fftw_execute(fplanN);
        
		// saving the calculated result
		for(int in1 = 0; in1 <nump; in1 ++)
            for(int in2 = 0; in2 <nump2; in2 ++)
                result[in1][in2] = complexv( FFTdatO[in1*nump2+in2][0], FFTdatO[in1*nump2+in2][1] );
	}
	else
	{
		cleanFFTSpace();
		executeN(result,source,in1,in2);
	}
}


interpolationF<complexv> toolsFFT::executeP(interpolationF<complexv>& source)
{
    // initial value is always assumed to be zero
    constants cst;
    double omd = source.GetXF();
    int omn = source.GetN();
    //cout<<omn<<" omn\n";
    omd = cst.pi2/omd;
    interpolationF<complexv> result(source);
    complexv* src = source.DirectAccessD();
    complexv* res = result.DirectAccessD();
    executeP(res,src, omn);
    result.UpdateAxis(0,omd);
    return result;
}
interpolationF<complexv> toolsFFT::executeN(interpolationF<complexv>& source)
{
    // initial value is always assumed to be zero
    constants cst;
    double omd = source.GetXF();
    int omn = source.GetN();
    //cout<<omn<<" omn\n";
    omd = cst.pi2/omd;
    interpolationF<complexv> result(source);
    complexv* src = source.DirectAccessD();
    complexv* res = result.DirectAccessD();
    executeN(res,src, omn);
    result.UpdateAxis(0,omd);
    return result;
}
void toolsFFT::SwapSides(interpolationF<complexv>& funct)
{
    // moves the dataset cyclicaly
    // update initial value
    int num = funct.GetN();
    complexv* res = new complexv[num];
    int nu2 = num/2;
    complexv* src = funct.DirectAccessD();

    for(int ind = 0; ind<nu2; ind++)
    {
        res[ind+nu2] = src[ind];
    }
    for(int ind = nu2; ind<num; ind++)
    {
        res[ind-nu2] = src[ind];
    }
    funct.DirectExchange(res);
    double max = funct.GetXF();
    double min = funct.GetXI();
    double st = (max-min)/num;
    if(min == 0)
            funct.UpdateAxis(-max/2,st);
    else
            funct.UpdateAxis(0,st);
    delete[] src;
}


// two dimensional sophisticated
interpolationF2d<complexv> toolsFFT::executeP(interpolationF2d<complexv>& source)
{

    // source initial values are assumed to be zero
    constants cst;

    double omdx = source.GetXF();
    double omdy = source.GetYF();
    int omnx = source.GetNx();
    int omny = source.GetNy();
    
    omdx = cst.pi2/omdx;
    omdy = cst.pi2/omdy;
    interpolationF2d<complexv> result(source);
    complexv** src = source.DirectAccessD();
    complexv** res = result.DirectAccessD();
    executeP(res,src, omnx,omny);
    result.UpdateAxis(0,omdx,0,omdy);
    return result;
}
interpolationF2d<complexv> toolsFFT::executeN(interpolationF2d<complexv>& source)
{

    // source initial values are assumed to be zero
    constants cst;

    double omdx = source.GetXF();
    double omdy = source.GetYF();
    int omnx = source.GetNx();
    int omny = source.GetNy();
    
    omdx = cst.pi2/omdx;
    omdy = cst.pi2/omdy;
    interpolationF2d<complexv> result(source);
    complexv** src = source.DirectAccessD();
    complexv** res = result.DirectAccessD();
    executeN(res,src, omnx,omny);
    result.UpdateAxis(0,omdx,0,omdy);
    return result;
}
void toolsFFT::SwapSides(interpolationF2d<complexv>& funct)
{
    // moves the dataset cyclicaly
    // update initial value
    int numx = funct.GetNx();
    int numy = funct.GetNy();
    int nu2x = numx/2;
    int nu2y = numy/2;
    
    storage<complexv> sres(2);
    sres.Allocate(numx,numy);
    complexv** res = sres.data2D;
    complexv** src = funct.DirectAccessD();

    // dimension 1
    for(int inx = 0; inx<numx; inx++)
    {
        for(int iny = 0; iny<nu2y; iny++)
        {
                res[inx][iny+nu2y] = src[inx][iny];
        }
        for(int iny = nu2y; iny<numy; iny++)
        {
                res[inx][iny-nu2y] = src[inx][iny];
        }
    }
    // dimension 2
    for(int iny = 0; iny<numy; iny++)
    {
        for(int inx = 0; inx<nu2x; inx++)
        {
                src[inx+nu2x][iny] = res[inx][iny];
        }
        for(int inx = nu2x; inx<numx; inx++)
        {
                src[inx-nu2x][iny] = res[inx][iny];
        }
    }
    
    sres.Delete();
    
    double maxx = funct.GetXF();
    double minx = funct.GetXI();
    double stx = (maxx-minx)/numx;
    double maxy = funct.GetYF();
    double miny = funct.GetYI();
    double sty = (maxy-miny)/numy;
    if(minx == 0 && miny == 0 )
            funct.UpdateAxis(-maxx/2,stx, -maxy/2,sty );
    else if(minx != 0 && miny != 0 )
            funct.UpdateAxis(0,stx, 0, sty);
    else if(minx == 0 && miny != 0 )
            funct.UpdateAxis(-maxx/2,stx, 0, sty);
    else
            funct.UpdateAxis(0,stx, -maxy/2, sty);
}
