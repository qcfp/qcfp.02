#include"interpolationF2d.hpp"



template<> 
double interpolationF2d<double>::Get(double x,double y)
{
    if(xstep<0 || ystep<0)
        return GetNegative(x,y);
        
    

        if(num_samplx == 0 || num_samply == 0)
        {
            cout<<"Error: interpolationF2d<valT> undefined\n";
            return 0.0;
        }

        else if(num_samplx == 1 && num_samply ==1)
        {
            return values->data2D[0][0];
        }

        
        
        else if(x<xmin)
        {
            // simple linear interpolation on y axis
            //return values->data2D[0][y];
            
            double* inX = new double[num_samply];
            double* inY = new double[num_samply];
            
            for(int ind = 0; ind<num_samply; ind ++)
            {
                inX[ind] = ymin + ystep*ind;
                inY[ind] = values->data2D[0][ind];
            }
            
            double ret, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inY,num_samply,y,ret,err);

            return ret;
        }
        
        else if(x>(xmax-xstep))
        {
            //return values->data1D[num_samplx-1][y];
           // simple linear interpolation on y axis
            
            double* inX = new double[num_samply];
            double* inY = new double[num_samply];
            
            for(int ind = 0; ind<num_samply; ind ++)
            {
                inX[ind] = ymin + ystep*ind;
                inY[ind] = values->data2D[num_samplx-1][ind];
            }
            
            double ret, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inY,num_samply,y,ret,err);

            return ret;

        }
    
        else if(y<ymin)
        {
            // simple linear interpolation on x axis
            //return values->data2D[x][0];
            
            double* inX = new double[num_samplx];
            double* inY = new double[num_samplx];
            
            for(int ind = 0; ind<num_samplx; ind ++)
            {
                inX[ind] = xmin + xstep*ind;
                inY[ind] = values->data2D[ind][0];
            }
            
            double ret, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inY,num_samplx,x,ret,err);

            return ret;

        }
        
        else if(y>(ymax-ystep))
        {
            // simple linear interpolation on x axis
            //return values->data2D[x][max];
            
            double* inX = new double[num_samplx];
            double* inY = new double[num_samplx];
            
            for(int ind = 0; ind<num_samplx; ind ++)
            {
                inX[ind] = xmin + xstep*ind;
                inY[ind] = values->data2D[ind][num_samplx-1];
            }
            
            double ret, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inY,num_samplx,x,ret,err);

            return ret;

        }
    
        
        // interpolation in 2D
        
        // getting index x
        double intervalx = xmax-xmin;
	double zx=(x-xmin)/intervalx*num_samplx;
	int iix= (int)floor(zx);
        
        // getting displacement from ii
	zx=zx-(double)iix;
	
        
        // getting index y
        double intervaly = ymax-ymin;
	double zy=(y-ymin)/intervaly*num_samply;
	int iiy= (int)floor(zy);
        
        // getting displacement from ii
	zy=zy-(double)iiy;
	
        
        
        // interpolating by a square
        double val4[4];
        val4[0]=values->data2D[iix][iiy];
        val4[1]=values->data2D[iix][iiy+1];
        val4[2]=values->data2D[iix+1][iiy+1];
        val4[3]=values->data2D[iix+1][iiy];
        
        toolsInterpolate2d intobj;
        return intobj.GetSquare0011(zx,zy,val4);
        
}


template<> 
complexv interpolationF2d<complexv>::Get(double x,double y)
{

    if(xstep<0 || ystep<0)
        return GetNegative(x,y);
        

    if(num_samplx == 0 || num_samply == 0)
        {
            cout<<"Error: interpolationF2d<valT> undefined\n";
            return 0.0;
        }

        else if(num_samplx == 1 && num_samply ==1)
        {
            return values->data2D[0][0];
        }

        
        
        else if(x<xmin)
        {
            // simple linear interpolation on y axis
            //return values->data2D[0][y];
            
            double* inX = new double[num_samply];
            double* inYR = new double[num_samply];
            double* inYI = new double[num_samply];
            
            for(int ind = 0; ind<num_samply; ind ++)
            {
                inX[ind] = ymin + ystep*ind;
                inYR[ind] = (values->data2D[0][ind]).real();
                inYI[ind] = (values->data2D[0][ind]).imag();
            }
            
            double retR, retI, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inYR,num_samply,y,retR,err);
            iobj.linint(inX,inYI,num_samply,y,retI,err);

            return complexv(retR,retI);
        }
        
        else if(x>=(xmax-xstep))
        {
            //return values->data1D[num_samplx-1][y];
           // simple linear interpolation on y axis
            
            double* inX = new double[num_samply];
            double* inYR = new double[num_samply];
            double* inYI = new double[num_samply];
            
            for(int ind = 0; ind<num_samply; ind ++)
            {
                inX[ind] = ymin + ystep*ind;
                inYR[ind] = (values->data2D[num_samplx-1][ind]).real();
                inYI[ind] = (values->data2D[num_samplx-1][ind]).imag();
            }
            
            double retR, retI, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inYR,num_samply,y,retR,err);
            iobj.linint(inX,inYI,num_samply,y,retI,err);

            return complexv(retR,retI);

        }
    
        else if(y<ymin)
        {
            // simple linear interpolation on x axis
            //return values->data2D[x][0];
            
            double* inX = new double[num_samplx];
            double* inYR = new double[num_samplx];
            double* inYI = new double[num_samplx];
            
            for(int ind = 0; ind<num_samplx; ind ++)
            {
                inX[ind] = xmin + xstep*ind;
                inYR[ind] = (values->data2D[ind][0]).real();
                inYI[ind] = (values->data2D[ind][0]).imag();
            }
            
            double retR, retI, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inYR,num_samplx,y,retR,err);
            iobj.linint(inX,inYI,num_samplx,y,retI,err);

            return complexv(retR,retI);


        }
        
        else if(y>=(ymax-ystep))
        {
            // simple linear interpolation on x axis
            //return values->data2D[x][max];
            
            double* inX = new double[num_samplx];
            double* inYR = new double[num_samplx];
            double* inYI = new double[num_samplx];
            
            for(int ind = 0; ind<num_samplx; ind ++)
            {
                inX[ind] = xmin + xstep*ind;
                inYR[ind] = (values->data2D[ind][num_samply-1]).real();
                inYI[ind] = (values->data2D[ind][num_samply-1]).imag();
            }
            
            double retR, retI, err;
            toolsInterpolate iobj;
            iobj.linint(inX,inYR,num_samplx,y,retR,err);
            iobj.linint(inX,inYI,num_samplx,y,retI,err);

            return complexv(retR,retI);

        }
    
        
        // interpolation in 2D
        
        // getting index x
        double intervalx = xmax-xmin;
	double zx=(x-xmin)/intervalx*num_samplx;
	int iix= (int)floor(zx);
        
        // getting displacement from ii
	zx=zx-(double)iix;
	
        
        // getting index y
        double intervaly = ymax-ymin;
	double zy=(y-ymin)/intervaly*num_samply;
	int iiy= (int)floor(zy);
        
        // getting displacement from ii
	zy=zy-(double)iiy;
	
        
        
        // interpolating by a square
        double val4[4];
        val4[0]=(values->data2D[iix][iiy]).real();
        val4[1]=(values->data2D[iix][iiy+1]).real();
        val4[2]=(values->data2D[iix+1][iiy+1]).real();
        val4[3]=(values->data2D[iix+1][iiy]).real();
        
        toolsInterpolate2d intobj;
        double retR,retI;
        retR = intobj.GetSquare0011(zx,zy,val4);
        
        val4[0]=(values->data2D[iix][iiy]).imag();
        val4[1]=(values->data2D[iix][iiy+1]).imag();
        val4[2]=(values->data2D[iix+1][iiy+1]).imag();
        val4[3]=(values->data2D[iix+1][iiy]).imag();
        
        retI = intobj.GetSquare0011(zx,zy,val4);

        return complexv(retR,retI);
}
