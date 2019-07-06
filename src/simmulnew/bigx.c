/************************************************************
*
*     COPYRIGHT (c) M.Bos 1998, for details see file copying
*
*     bigx.c
*     calculation of double integral
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

double bigx(double p, double q, double d1, double d2, double mu1j,
            double mu2j, double mubj, double db)
{
   double vd1d2, vd10, v0d2, v00;
   double t1, t2, t3, result ;    
  
   if  ( db==0.0)
    {
      t1 = -1.0/(p*mu1j+q*mu2j);
      t2 = (mu2j/p)*log(fabs(1.0 + p/mu2j));
      t3 = (mu1j/q)*log(fabs(1.0 - q/mu1j)) ;
      v00 = t1*( t2 + t3);

      /* debugging  
  printf("in bigx with db=0 t1 %g t2 %g t3 %g v00 %g\n", t1, t2, t3, v00); /* */

    }
   else
    {
      v00 = bigv(0.0, 0.0, mu1j, mu2j, mubj, db, p, q);

      /* debugging  
      printf("In bigx db <> 0.0 v00 is %g \n", v00); /* */

    }

   vd1d2 = bigv(d1,d2, mu1j, mu2j, mubj, db, p, q);
   vd10 = bigv(d1, 0.0, mu1j, mu2j, mubj, db, p, q);
   v0d2 = bigv(0.0, d2, mu1j, mu2j, mubj, db, p, q);

   /* debugging  
   printf(" In bigx vd1d2 %g vd10 %g v0d2 %g v00 %g \n",
           vd1d2, vd10, v0d2, v00); /* */

   result=(vd1d2 - vd10 - v0d2 + v00);

   /* debugging 
   printf("bigx is %g \n", result); /* */

   return(result);
}


double bigv(double d1, double d2, double mu1j, double mu2j, double mubj,
            double db, double p, double q)
{
  double t1, t2, t3, t4, t5, t6 ;

  /* debugging 
  printf("In bigv p %g q %g mu1j %g mu2j %g mubj %g db %g \n",
          p, q, mu1j, mu2j, mubj, db);
  printf("d1 is %g d2 is %g \n", d1, d2);  /* */

  t1 = exp( (q-mu1j)*d1 -(p+mu2j)*d2 - mubj*db) ;
  
  t2 = mu2j/( p*(p*mu1j + q*mu2j)) ;
  t3 = (1.0 + p/mu2j)*(mu1j*d1 + mubj*db + mu2j*d2);

  /* debugging  

  printf("t3 voor expon integral calc %g \n", t3);  /* */

  t3 = dfunc(t3);
  t4 = mu1j/(q*(p*mu1j+q*mu2j));
  t5 = (1.0 - q/mu1j)*(mu1j*d1 + mubj*db + mu2j*d2) ;

  /* debugging 
  printf("In bigv t5 voor exp int %g \n", t5); /* */

  t5 = dfunc(t5);
  t6 = mu1j*d1 + mubj*db + mu2j*d2 ;

  /* debugging 
  printf("t6 voor expon integral calc %g \n", t6);  /* */

  t6 = dfunc(t6);
  t6= (-1.0/(p*q))*t6 ;

  /* debugging 
  printf("In bigv === t1 %g t2 %g t3 %g t4 %g t5 %g t6 %g \n",
           t1, t2, t3, t4, t5, t6);  /* */

  return t1*(t2*t3 + t4*t5 + t6);
}

