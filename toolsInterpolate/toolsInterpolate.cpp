#include"toolsInterpolate.hpp"
#include"cmath"
#include<iostream>

toolsInterpolate::toolsInterpolate()
{
    tint = 0;
}
toolsInterpolate::~toolsInterpolate()
{
    interclean();
}
void toolsInterpolate::interclean()
{
    if(tint != 0)
        delete[] tdata;
    tint = 0;
}


// rational function interpolation
int toolsInterpolate::ratint(double* xa, double* ya, int n, double x, double& y, double& dy)
{
	int m,i,ns=0;
	double w,t,hh,h,dd;
	
	double* c=new double[n];
	double* d=new double[n];
	hh=fabs(x-xa[0]);
	for(i=0;i<n;i++) 
    {
		h=fabs(x-xa[i]);
		if (h == 0.0 ) {
			dy=0.0;
            delete[] c;
            delete[] d;
			y = ya[i];
            return 0;
		} 
        else if (h < hh) 
        {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]*1.000000000000001;
	}
	y=ya[ns--];
	for(m=1; m<n;m++) 
    {
		for(i=0; i<(n-m) ;i++) 
        {
			w=c[i+1]-d[i];
			h=xa[i+m]-x;
			t= (xa[i]-x)*d[i]/h;
			dd=t-c[i+1];
			if (dd == 0.0 ) 
            {
                delete[] c;
                delete[] d;
                y = 0;
                dy = ya[0];
                return 1;
                //"Error: this is a pole";
            }
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		if(2*(ns+1) < (n-m))
			dy= c[ns+1];
		else dy= d[ns--];
		y += dy;
	}
    delete[] c;
    delete[] d;
    return 0;
}


// polynomial interpolation
int toolsInterpolate::polint(double *inX, double *inY, int n, double x,  double& ret, double& err)
{
    int i,j,k;
    double t1,tt1;
    
	double* c=new double[n];
	double* d=new double[n];
    
    err=0;
    
    t1=fabs(x-inX[0]);
    k=0;
    //initial values
    for(i=0;i<n;i++)
    {
        c[i]=inY[i];
        d[i]=inY[i];
        
        if(x==inX[i]) 
        {
            ret = inY[i];
            delete[] c;
            delete[] d;
            err = 0.0;
            return 0;
        }
        
        //now I will search for the nearest result
        if((tt1=fabs(x-inX[i]))<t1)
        {t1=tt1;k=i;}
        
    }
    t1=c[k];
    for(i=1;i<n;i++)
    {
        
        for(j=0;j<n-i;j++)
        {
            tt1 = (c[j+1]-d[j])/(inX[j]-inX[j+i]);
            c[j]=(inX[j]-x)*tt1;
            d[j]=(inX[j+i]-x)*tt1;
        }
        
        if(2*k<(n-i)){t1+=c[k];err=c[k];}
        else {t1+=d[k-1];k--;err=d[k];}
    }
        
    ret = t1;
    
    delete[] c;
    delete[] d;
    return 0;
}

// Barycentric interpolation (numerical recipies)
int toolsInterpolate::barint(double *inX, double *inY, int n, double x,  double& ret, double& err)
// this default setup uses d=n-1 interpolation
{
    int retval;
    retval = barinit(inX,n,n-1);
    if(retval == 0)
        retval = barretval(inX, inY, x, ret);
    interclean();
    err = 0; // not used
    return retval;

}

int toolsInterpolate::barinit(double *inX, int n, int d)
// this function calculates the weights
// for the d order of Barycentric interpolation
{
    
    if (n<=d) 
    {
        std::cout<<"Barycentric interpolation order too large for number of points";
        return 1;
    }
    
    interclean();
    tdata =  new double[n];
    tint = n;
    
    for (int k=0;k<n;k++) 
    {
        int imin= (k-d<0? 0: k-d);
        int imax = (k >= n-d? n-d-1: k);
        double temp = (imin & 1? -1.0: 1.0);
        double sum=0.0;
        for(int i=imin;i<=imax;i++)
        { 
            int jmax=(i+d<n-1 ? i+d: n-1);
            double term=1.0;
            for (int j=i;j<=jmax;j++) 
            {
                if (j==k) continue;
                term *= (inX[k]-inX[j]);
            }
            term=temp/term;
            temp=-temp;
            sum += term;
        }
        tdata[k]=sum; 
    }
    return 0;
}
int toolsInterpolate::barretval(double *inX, double *inY, double x, double& y)
{
    int n = tint;
    double num=0,den=0;
    for (int i=0;i<n;i++) 
    {
        double h=x-inX[i];
        if (h == 0.0) {
            y = inY[i];
            return 0;
        } 
        else 
        {
            double temp=tdata[i]/h;
            num += temp*inY[i];
            den += temp;
        } 
    }
    y = num/den;
    return 0;
}


// simple linear interpolation
int toolsInterpolate::linint(
                             double *inX,
                             double *inY,
                             int n,
                             double x,
                             double& ret,
                             double& err)
{
    err = 0;
    
    int ix;
    int cse;
    // looking for index corresponding to x
    for(ix=0; ix<n; ix ++)
    {
        if(x<inX[ix])
        {
            cse = 1;
            break;
        }
        else if(x>inX[n-ix-1])
        {
            cse = 0;
            break;
        }
    }
    
    if(cse)
    {
        if(ix == 0)
        {
            ret = inY[0];
            return 0;
        }
        else
            ix --;
    }
    else // cse == 0
    {
        
        if(ix == 0)
        {
            ret = inY[n-1];
            return 0;
        }
        else
            ix = n-ix-1;
    }
    
    ret= (inY[ix+1]-inY[ix])/(inX[ix+1]-inX[ix])*(x-inX[ix])+inY[ix];
    
    return 0;
    
}
