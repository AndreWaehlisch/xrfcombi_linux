/************************************************************
*
*     COPYRIGHT (c) M.Bos 1998, for details see file copying
*
*     mupdown.c
*     calculation of Mup and Mdown terms
*     cf D.K.G. de Boer X-Ray Spectrom. 19 (1990) 145
*
*     M.Bos
*     sept 1998
*
************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

double mdown(double mu1lprim, double mu1iprim, double mublprim,
             double mu2lprim, double mu1j, double mu2j, double mubj,
             double d1, double d2 , double db)
{
   double t1, t2, t3 ;

   t1 = (mu1lprim + mu1iprim)/(1.0 - exp(-(mu1lprim + mu1iprim)*d1)) ;
   if (db == 0.0)
     t2 = 1.0 ;
   else
      t2 = exp(mublprim*db) ;
   t3 = bigx(mu1iprim, mu2lprim, d2, d1, mu1j, mu2j, mubj, db);
   
   return (t1*t2*t3*exp(100));
}

double mup(double mu1lprim, double mu1iprim, double mublprim,
             double mu2lprim, double mu1j, double mu2j, double mubj, 
             double d1, double d2 , double db)
{
   double t1, t2, t3 ;

   t1 = (mu1lprim + mu1iprim)/(exp((mu1lprim + mu1iprim)*d1) - 1.0) ;
   if (db == 0.0)
     t2=1.0 ;
   else
      t2 = exp(-mublprim*db) ;

   t3 = bigx(mu2lprim, mu1iprim, d1, d2, mu1j, mu2j, mubj, db);

   /* debugging
   printf("In mup t1 %g t2 %g t3 %g \n", t1, t2, t3); /* */
   
   return (t1*t2*t3*exp(100));
}
