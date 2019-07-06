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
    
  double ps2[9]={ 1.737331760720576030932E-8,
                  9.999989642347613068437E-1,
                  1.487967702840464066613E1,
                  7.633628843705946890896E1,
                  1.698106763764238382705E2,
                  1.700632978311516129328E2,
                  7.246689782858597021199E1,
                  1.107326627786831743809E1,
                  3.828573121022477169108E-1};

  double qs2[9]={ 1.000000000000000000E00,
                  1.587964570758947927903E1,
                  9.021658450529372642314E1,
                  2.342573504717625153053E2,
                  2.953136335677908517423E2,
                  1.775728186717289799677E2,
                  4.662179610356861756812E1,
                  4.344836335509282083360E0,
                  8.258160008564488034698E-2};
                    
  double ps3[10]={ -9.9999999999999999087819E-1,
                   -5.2199632588522572481039E1,
                   -1.0611777263550331766871E3,
                   -1.0816852399095915622498E4,
                   -5.9346841538837119172356E4,
                   -1.7503273087497081314708E5,
                   -2.6181454937205639647381E5,
                   -1.7283375773777593926828E5,
                   -3.5846198743996904308695E4,
                   -1.3276881505637444622987E2};

  double qs3[10]= { 1.0000000000000000000000E0,
                    5.4199632588522559414924E1,
                    1.1635769915320848035459E3,
                    1.2842808586627297365998E4,
                    7.9231787945279043698718E4,
                    2.7858134710520842139357E5,
                    5.4616842050691155735758E5,
                    5.5903756210022864003380E5,
                    2.5989762083608489777411E5,
                    3.9147856245556345627078E4};



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
     for (i=0; i <=8 ; i++)
     {
       numerator+= ps2[i]*pow(x,-i);
       denominator += qs2[i]*pow(x,-i);
     }
     result *=numerator/denominator;
     return(result);
   }
   
   numerator=0.0;
   denominator=0.0;
   for (i=0; i <=9 ; i++)
   {
     numerator+= ps3[i]*pow(x,-i);
     denominator += qs3[i]*pow(x,-i);
   }
   result = (exp(-x)/x)*(1.0+  numerator/(x*denominator));
   return(result);
}

double dfunc(double x )
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
    
  double ps2[9]={ 1.737331760720576030932E-8,
                  9.999989642347613068437E-1,
                  1.487967702840464066613E1,
                  7.633628843705946890896E1,
                  1.698106763764238382705E2,
                  1.700632978311516129328E2,
                  7.246689782858597021199E1,
                  1.107326627786831743809E1,
                  3.828573121022477169108E-1};

  double qs2[9]={ 1.000000000000000000E00,
                  1.587964570758947927903E1,
                  9.021658450529372642314E1,
                  2.342573504717625153053E2,
                  2.953136335677908517423E2,
                  1.775728186717289799677E2,
                  4.662179610356861756812E1,
                  4.344836335509282083360E0,
                  8.258160008564488034698E-2};
                    
  double ps3[10]={ -9.9999999999999999087819E-1,
                   -5.2199632588522572481039E1,
                   -1.0611777263550331766871E3,
                   -1.0816852399095915622498E4,
                   -5.9346841538837119172356E4,
                   -1.7503273087497081314708E5,
                   -2.6181454937205639647381E5,
                   -1.7283375773777593926828E5,
                   -3.5846198743996904308695E4,
                   -1.3276881505637444622987E2};

  double qs3[10]= { 1.0000000000000000000000E0,
                    5.4199632588522559414924E1,
                    1.1635769915320848035459E3,
                    1.2842808586627297365998E4,
                    7.9231787945279043698718E4,
                    2.7858134710520842139357E5,
                    5.4616842050691155735758E5,
                    5.5903756210022864003380E5,
                    2.5989762083608489777411E5,
                    3.9147856245556345627078E4};



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
        result *= exp(x);
      return(result);
    }
   if ( x <= 4.0)
   {
     result= 1 ;
     numerator=0.0;
     denominator=0.0;
     for (i=0; i <=8 ; i++)
     {
       numerator+= ps2[i]*pow(x,-i);
       denominator += qs2[i]*pow(x,-i);
     }
     result *=numerator/denominator;
     return(result);
   }
   
   numerator=0.0;
   denominator=0.0;
   for (i=0; i <=9 ; i++)
   {
     numerator+= ps3[i]*pow(x,-i);
     denominator += qs3[i]*pow(x,-i);
   }
   result = (1.0/x)*(1.0+  numerator/(x*denominator));
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
/*
double dfunc(double x)
{
      return(exp(x)*E1z(x));
}
*/ 

 

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
