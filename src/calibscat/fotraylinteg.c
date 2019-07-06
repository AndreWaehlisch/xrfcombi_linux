/*
*    file fotraylinteg.c
*
*    functions to calculate integrand
*    for double integral for crosssection
*    rayleigh scattering
*
*    uses a numer of global parameters,
*    i.e. 
*    double EEi
*    double muci
*    double muclabda
*    layer *ptrlaag
*    double sinpsi1
*    double sinpsi2
*    double cospsi2
*
*    and two dependent variables alfa and fi
*   

*    Dec 12th 2003 changed integr_func1 and 2
*    removed factor sin(a) resp -sin(alfa) 
*    cf paper tirao et al XRS 32, 2003, p 13-24
*    without sin(alfa) in numerator
*
*    Dec 13th sin(alfa) moved back in
*    monochromatic excitation fits better
*    with sin(alfa) in.
*
*    Changed to use Tirao et al eqns
*
*    code by M.Bos
*
*    Dec 19th 2003
*
*/

#include <stdio.h>
#include <math.h>
#include "xrfluor.h"

extern double EEi ;
extern double muci ;
extern double muclabda;
extern layer *ptrlaag;
extern double sinpsi1;
extern double sinpsi2;
extern double cospsi1;
extern double cospsi2;
extern double store_thickness ;
extern elem_list lijst_els;
extern int f2cnt;


double fotrayl_func1(double alfa, double fi)
{
  double theta, costheta, cross, result, thetadegr;
  double part1, part2, part3 ;
  double thickness ;
  /* calc theta for case t1 < t2 */
  
  costheta= -sinpsi2*cos(alfa)-sin(alfa)*cos(fi)*cospsi2;
  theta= acos(costheta);
  thetadegr= (theta/PI)*180.0;
  /* debugging 
  printf("in f1 theta is %g \n",theta); /* */
  /* debugging 
  printf("in f1 alfa %g fi %g theta %g sinpsi2 %g cospsi2 %g \n", 
         alfa,fi,theta, sinpsi2, cospsi2); /* */
    thickness = store_thickness ;

    cross=raycrossmix(EEi, ptrlaag, &lijst_els, thetadegr); /* */
/*	 printf("in f1 cross is %g bij thetadegr %g\n", cross, thetadegr); /*  */

   part1= 1/(muci+cos(alfa)*muci/sinpsi2);
 /*  printf("f1 part1 %e \n",part1);  /*  */
   
   part2 = exp(-(muclabda/sinpsi1+muci/sinpsi2)*thickness);
   part2 /= (cos(alfa)*muclabda/sinpsi1 -muci) ;
 /*  printf("f1 part2 %e\n",part2);  /*  */

   part3= -cos(alfa)*exp(-(muci/cos(alfa)+muci/sinpsi2)*thickness);
   part3 *= (muclabda/sinpsi1+muci/sinpsi2) ;  /* is in denom pil */
   part3 /= (muci+cos(alfa)*muci/sinpsi2) ;
   part3 /= (cos(alfa)*muclabda/sinpsi1-muci);

/*   printf("f1 part3 %e\n",part3);  /*  */

   result= part1 + part2 + part3 ;   /*  */

   result *= cross*sin(alfa) ;     /* sin a added pa0,bo */
 
  /* printf("in f1 result is %g alfa %g fi %g thetadegr %g\n", result, alfa, fi, thetadegr); /* */
  return(result);
}


double fotrayl_func2(double alfa, double fi)
{
  double theta, costheta, cross, result, thetadegr;
  double part1, part2, part3 ;
  double thickness ;

  /* calc theta for case t1 < t2 */
 
 /* if (fabs(alfa - PI/2)< 1.0e-10)
	   return(0.0);   */

   
  costheta= -sinpsi2*cos(alfa)+sin(alfa)*cos(fi)*cospsi1;
  theta= acos(costheta);
   thetadegr= (theta/PI)*180.0 ;
  /* debugging 
  printf("in f2 alfa %g fi %g theta %g sinpsi2 %g cospsi1 %g \n", 
         alfa,fi,theta, sinpsi2, cospsi1); /* */
  thickness = store_thickness;   
  cross=raycrossmix(EEi, ptrlaag, &lijst_els, thetadegr); /*  */
 /*  printf("in f2  alfa is %g  cos %g tan %g , fi %g \n",alfa,  cos(alfa), tan(alfa), fi); /* */
  /* printf("in f2 cross is %g bij theta %g\n,", cross, thetadegr); /* */

   part1= 1.0 / ( muci-cos(alfa)*muclabda/sinpsi1);

  /* printf("f2 part1  %g \n", part1);  /* */


   part2= exp(-(muclabda/sinpsi1+muci/sinpsi2)*thickness);
 /* printf("f2 part2 step 1 %g\n",part2);  /*  */
   part2 /= (-cos(alfa)*muci/sinpsi2 -muci);
 /*  printf("f2 part2 step2 %g alfa %g muci %g sinpsi2 % g \n",
	   part2, alfa, muci, sinpsi2); /* */
   
   part3 =   cos(alfa)*exp(-(muclabda/sinpsi1 - muci/cos(alfa))*thickness);
  /* printf("f2 part3 step 1 %g \n",part3);  /* */
   part3 *=  (muclabda/sinpsi1+muci/sinpsi2);
   part3 /= (muci+cos(alfa)*muci/sinpsi2);
 /* printf("f2 part3 step 3 %g alfa %g muci %g sinpsi2 %g fi %g\n",
	       part3, alfa, muci, sinpsi2, fi); /*  */
   part3 /= (cos(alfa)*muclabda/sinpsi1 - muci);
 /* printf("f2 part3 step 4 %g\n",part3);  /* */

   
   result= sin(alfa)*cross*(part1 +  part2 + part3) ;
 
   /* printf("in f2 result is %g\n", result); /* */
 /* printf(" in f2 muci is %g muclabda is %g\n", muci, muclabda); /* */ 
  f2cnt++ ;  /* registratie functiecalls */
  return(result);
}

