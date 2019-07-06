/****************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
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
*
*  exponential integral acc. W.J. Cody, H.C. Thacher Jr
*  Math. Comput. 22 (1968) 641
*
*  M.Bos
*  sept 1997
*
******************************************************/
#include <stdio.h>
#include <math.h>

double E1z(double );
double Eiz(double, int );

double dfunc(double);
double lzero(double, double, double, double);
double ltot(double, double, double, double);
double linfty(double, double, double);

double E1z(double x )
{
    int i ;
    double result, numerator, denominator;

    double ps1[7]={ -1.4815102102575750838086E05,
                    1.5026059476436982420737E05,
                    8.9904972007457256553251E04,
                    1.5924175980637303639884E04,
                    2.1500672908092918123209E03,
                    1.1669552669734461083368E02,
                    5.0196785185439843791020E00} ;
   double qs1[7]={ 2.5666493484897117319268E05,
                   1.8434070063353677359298E05,
                   5.2440529172056355429883E04,
                   8.1258035174768735759866E03,
                   7.5043163907103936624165E02,
                   4.0205465640027706061433E01,
                   1.0000000000000000000000E00} ; 
    
  double ps2[8]={ 8.677459548384437437E-8,
                  9.999955193013903006E-1,
                  1.184831055549458443E01,
                  4.559306442533898233E01,
                  6.992794512910030229E01,
                  4.252020347688407791E01,
                  8.836718088038439386E00,
                  4.013776649406647203E-1};

  double qs2[8]={ 1.000000000000000000E00,
                  1.284819353791566499E01,
                  5.644335695618031986E01,
                  1.066451837699138825E02,
                  8.973110971252898022E01,
                  3.149718491704407502E01,
                  3.795590037621222428E00,
                  9.088045691888692188E-2};
   
  double ps3[8]={ -9.9999999999997341443E-1,
                  -3.4406199500668489491E01,
                  -4.2753267120198853941E02,
                  -2.3960194324749054028E03,
                  -6.1688521005547635088E03,
                  -6.5760969874802117925E03,
                  -2.1060773714263328896E03,
                  -1.4899084997294816902E01};

  double qs3[8]= { 1.0000000000000000000E00,
                   3.6406199500645980400E01,
                   4.9434507020990364527E02,
                   3.1902723748954330383E03,
                   1.0337075308584097698E04,
                   1.6324145355778350289E04,
                   1.1149775287109662000E04,
                   2.3781389910216022120E03};

 

    if ( x <= 1.0)
    {
      result=  -log(fabs(x)) ;
      numerator=0.0 ;
      denominator=0.0;
      for (i=0; i <= 6 ; i++)
      {
        numerator += ps1[i]*pow(x,i);
        denominator += qs1[i]*pow(x,i);
      }
        result += numerator/denominator;
      return(result);
    }
   if ( x <= 4.0)
   {
     result= exp(-x) ;
     numerator=0.0;
     denominator=0.0;
     for (i=0; i <=7 ; i++)
     {
       numerator+= ps2[i]*pow(x,-i);
       denominator += qs2[i]*pow(x,-i);
     }
     result *=numerator/denominator;
     return(result);
   }
   
   numerator=0.0;
   denominator=0.0;
   for (i=0; i <=7 ; i++)
   {
     numerator+= ps3[i]*pow(x,-i);
     denominator += qs3[i]*pow(x,-i);
   }
   result = (exp(-x)/x)*(1.0+  numerator/(x*denominator));
   return(result);
}

/*
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
          printf("som bij %d is %e \n", i, sum); 
    }
return(sum);
}
 */

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
      return(exp(x)*E1z(x));
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
