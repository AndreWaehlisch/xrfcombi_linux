/*
*	selfattilabda.c
*
*	calculates attenuation of primary fluorescence
*	due to thickness of the layer itself
*
*       D.K.G de Boer et.al. Adv. in X-ray Anal. 33 (1990) 237
*
*	M.Bos
*	august 1997
*/
#include <stdio.h>
#include <math.h>
#include "xrfluor.h"

extern double psi1;
extern double psi2;
 

double selfattn(double mu1lambda, double mu1i, double d)
{
   double macht;
   macht= -(mu1lambda/sin(psi1) + mu1i/sin(psi2))*d;
   return (1.0- exp(macht));
}
 
