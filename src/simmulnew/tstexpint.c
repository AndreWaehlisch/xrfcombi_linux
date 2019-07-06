/********************
*
* tstexpint.c
*
*   testroutine to drive calc exponential integral
*
*   M.Bos
*
*   sept 1998
*
*******************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double E1z(double);
double dfunc(double);

int  main(int argc, char *argv[])
{
   double argument ;

   argument= atof(argv[1]);

   printf("exp int van %s is %g dfunc is %g \n", argv[1], E1z(argument),
          dfunc(argument) );
   return 0 ;

}




