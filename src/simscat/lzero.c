/****************************************************
*
*  COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*  lzero.c
*
*  routine to calculate L0 factor needed
*  in calculation of secundary fluorescence
*
*  cf D.K.G. de Boer and P.N. Brouwer
*  Advances in X-ray Analysis
*  33 (1990) 237
*
*  M.Bos
*  sept 1997
*
******************************************************/
#include <stdio.h>
#include <math.h>

double E1z(double, int );
double Eiz(double, int );

double dfunc(double);
double lzero(double, double, double, double);
double ltot(double, double, double, double);
double linfty(double, double, double);

double E1z(double x, int n)
{
    int i ;
    double result, sigmaterm;
    result= -0.5772156649 - log(fabs(x));
    sigmaterm= -x;
    result -= sigmaterm;
    for (i=2; i <= n ; i++)
    {
        sigmaterm= (sigmaterm/i)*(-x);
        result -= sigmaterm/i;
    }
   return(result);
}

double reeks(double x)
{
    int i;    
    double sum, term, hulpterm;
    sum=0.0 ;
    term = x;

    for (i=2; i < 100 ; i++)
    {
          term=term*x;
          term=term/i;
          hulpterm=term/i;
          sum += hulpterm;
/*          printf("som bij %d is %e \n", i, sum); */
    }
return( sum);}
double Eiz(double x, int n)
{
    int i ;
    double result, sigmaterm;
    result= 0.5772156649 + log(fabs(x));
    sigmaterm= x;
    result += sigmaterm;
    for (i=2; i <= n ; i++)
    {
        sigmaterm= (sigmaterm/i)*x;
        result += sigmaterm/i;
    }
   return(result);
}

double dfunc(double x)
{
      return(exp(x)*E1z(x, 1000));
}
    

double lzero(double mu1, double mu2, double muj, double d)
{
  double mu1_2_d, expmu12d ;
  double t1, t2, t3, t4, t5, t6, t7;
  double result, mudt3, mudt4, mudt5  ;
  double mudlimit=40.0;
  mu1_2_d= (mu1 +mu2)*d ;
  expmu12d= exp(-mu1_2_d);
  mudt3= muj*d;
  if (mudt3 > mudlimit)
    {
    return(linfty(mu1,mu2,muj));
    }
  mudt4=fabs(muj-mu2)*d;
  if  (mudt4 > mudlimit)
   { 
    return(linfty(mu1,mu2,muj)); 
   }
  mudt5= (mu1+muj)*d;
  if (mudt5 > mudlimit)
    {
      return(linfty(mu1,mu2,muj));
    }
  t1= (mu1+mu2)*muj/(1-expmu12d);
  t2=exp(-(mu1+muj)*d);
  t3=-dfunc(mudt3)/(mu1*mu2);
  t4=dfunc((muj-mu2)*d)/(mu2*(mu1+mu2));
  t5=dfunc(mudt5)/(mu1*(mu1+mu2));
  t6=log(fabs(1.0+mu1/muj))/(mu1*(mu1+mu2));
  t7=log(fabs(1.0-mu2/muj))*expmu12d/(mu2*(mu1+mu2));
  result=t1*(t2*(t3+t4+t5)+t6+t7);
  return(result);
}

double ltot(double mu1, double mu2, double muj, double d)
{
   return(lzero(mu1,mu2,muj,d)+lzero(mu2,mu1,muj,d));
}

double linfty(double mu1, double mu2, double muj)
{
   double result; 
   result=(muj/mu1)*log(1.0+mu1/muj);
   return(result);
}
