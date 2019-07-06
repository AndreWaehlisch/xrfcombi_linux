/*
*   file errorf.c
*   M.Bos
*   calculates errorfunction complement
*   based on Chebyshev fitting
*   cf Numerical recipies in C p 221
*
*   March 23 2001
*
*/
#include <stdio.h>
#include <math.h>

#define LIMIT 5.0
#define AANTAL 100

double erfcc(double x);

int main(void)
{
   double x=0, step ;
   int i ;
   step = LIMIT/((double) AANTAL);

   for (i=0; i < AANTAL ; i++)
   {
      x = i*step ;
      printf("%lf %lg\n",x,erfcc(x));
   }
   return 0;
}

double erfcc(double x)
{
  double t,z,ans ;
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196 + t*(0.09678418+
      t*(-0.18628806 + t*(0.27886807 +t*(-1.13520398+t*(1.48851587+
      t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans ;
}

