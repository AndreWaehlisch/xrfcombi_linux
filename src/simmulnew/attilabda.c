/*
*	attilabda.c
*
*	calculates attenuation of primary fluorescence
*	due to layers on top  
*
*       D.K.G de Boer et.al. Adv. in X-ray Anal. 33 (1990) 237
*
*	M.Bos
*	august 1997
*/
#include <stdio.h>
#include <math.h>
#include "xrfluor.h"

 

double attn(double mualambda, double psi1, double muai,
                double psi2, double da)
{
   double macht;
   macht= -(mualambda/sin(psi1) + muai/sin(psi2))*da;
   return ( exp(macht));
}
 
