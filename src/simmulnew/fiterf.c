/*
*    file firerf.c
*    matches errorfunction factor 1/(2*sqrt(Dt))
*    on exponent found on fiiting XRF data
*
*    M.Bos
*    March 23  2001
*/
#include <stdio.h>
#include <math.h>

#define STEPS 1000
#define CHECKPTS 100 

double erfcc(double);

int main(void)
{
    double expfactor, erfactor, tmp, low, high, stepsize ;
    int i, jj ;
    double ssq ;
    double x;

   printf("Give low and high boundary of erfactor ");
    scanf("%lf %lf", &low, &high);
    printf("Give exponential decay factor  ");
    scanf("%lf", &expfactor);
    stepsize = (high-low)/STEPS ;  

    for (i=0 ; i < STEPS ; i++)
    {
      erfactor = low  + i*stepsize;
      ssq = 0.0 ;
      for ( jj = 0 ; jj < CHECKPTS ; jj++)
      {
         x= jj*(2.5e-3)/CHECKPTS;
         ssq += (exp(-expfactor*x) - erfcc(erfactor*x))*
                (exp(-expfactor*x) - erfcc(erfactor*x));
      }  
      printf("erfactor %lg  ssq %lg\n", erfactor, ssq); 
    }   
  return 0 ;
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

