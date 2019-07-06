/*
*    file raylfotinteg.c
*
*    functions to calculate integrand
*    for double integral for crosssection
*    rayleigh scattering
*
*    uses a number of global parameters,
*    i.e. 
*    double EEi
*    double muci
*    double muclabda
*    layer *ptrlaag
*    double sinpsi1
*    double sinpsi2
*    double cospsi2
*    double store_labda
*    double store_labdai
*
*    and two dependent variables alfa and fi
*   
*    now based on equations of Tirao et al XRS 32 (2003)13
*    
*
*    code by M.Bos
*
*    Jan 16th 2004
*
*    now uses cos(alfa) as integrand instead of alfa
*
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
extern elem_list lijst_els;
extern int f2cnt;
extern double store_labda ;
extern emudata abstabel ;
extern int store_z ;
extern double store_energy;
extern double store_thickness ;

double raylfot_func1(double cosa, double fi)
{
  double theta, costheta, cross, result, thetadegr;
  double  energy, part1, part2, part3, thickness, alfa ;
  
 /* printf("in func1 store_energy is %g \n", store_energy); /*  */
  /* calc theta for case t1 < t2 */
  alfa=acos(cosa) ;
/*  printf("raylf1 alfa is %g fi is %g cospis1 %g sinpsi1 %g\n", alfa, fi, cospsi1, sinpsi1);  /*  */
  costheta= sin(alfa)*cos(fi)*cospsi1+ cosa*sinpsi1;
/*  printf("costheta is %g \n",costheta); /* */
  theta= acos(costheta);
  energy =  store_energy;
  thickness = store_thickness ;
  
/*  printf("theta is %f \n",theta); /* */
  thetadegr= (theta/PI)*180.0;
  /* debugging */
 /* printf("in f1 theta is %g \n",theta); /* */
  /* debugging */
 /* printf("in f1 alfa %g fi %g costhe %g theta %g sinpsi2 %g cospsi2 %g PI %g\n", 
         alfa,fi, costheta,theta, sinpsi2, cospsi2, PI); /* */
    cross=raycrossmix(store_energy, ptrlaag, &lijst_els, thetadegr); /* */
  /* printf("in f1 cross is %g bij thetadegr %g\n", cross, thetadegr); /*  */
   part1 = 1/(muclabda + cosa*muci/sinpsi2) ;
/*   printf(" f1 part 1 %g \n",part1);  /* */
   part2 = exp(-(muclabda/sinpsi1 + muci/sinpsi2)*thickness);
   part2 /= (cosa*muclabda/sinpsi1 - muclabda);
/*   printf("f1 part 2 %g cosa %g muclabda %g muci %g \n", 
	       part2, cosa, muclabda, muci);  /*  */
   part3 = -cosa*exp(-(muclabda/cosa+muci/sinpsi2)*thickness);
   part3 *= (muclabda/sinpsi1 + muci/sinpsi2);
   part3 /= (muclabda + cosa*muci/sinpsi2)*(cosa*muclabda/sinpsi1 - muclabda);
/*   printf("f1 part3 %g \n", part3); /*  */

    result=(cross)*(part1+part2+part3) ;
/*	printf("in f1 result is %g at alfa %g and fi %g thetadegr %g\n", result, alfa, fi, thetadegr);  /* */
  return(result);
}


double raylfot_func2(double cosa, double fi)
{
  double theta, costheta, cross, result, thetadegr;
  double energy, part1, part2, part3, thickness, alfa ;
  alfa=acos(cosa);
  costheta= -sin(alfa)*cos(fi)*cospsi1+ cosa*sinpsi1;
  theta= acos(costheta);
  thetadegr= (theta/PI)*180.0 ;
  energy= store_energy;
  thickness = store_thickness ;
  /* debugging 
  printf("in f2 alfa %g fi %g theta %g sinpsi2 %g cospsi1 %g \n", 
         alfa,fi,theta, sinpsi2, cospsi1); /* */
  /* cross=0.1234 ; /* ter test cross opzoeken uitgeschakeld */
   cross=raycrossmix(store_energy, ptrlaag, &lijst_els, thetadegr); /*  */
 /*  printf("in f2  alfa is %g  cos %g tan %g , fi %g \n",alfa,  cos(alfa), tan(alfa), fi); /* */
  /* printf("in f2 cross is %g bij theta %g\n,", cross, thetadegr); /* */
  part1 = 1/(muclabda -cosa*muclabda/sinpsi1); 
/*  printf("f2 part1 %g \n",part1);  /*  */
  part2 = exp(-(muclabda/sinpsi1 + muci/sinpsi2)*thickness);
  part2 /= (-cosa*muci/sinpsi2 - muclabda);
 /* printf("f2 part 2 %g \n",part2);   /*  */
  
 /* printf("f2 part3a %g  cosa %g muclabda %g sinpsi1 %g d %g \n",
	   part3, cosa, muclabda, sinpsi1, thickness);  /*  */
  part3 = (muclabda/sinpsi1 + muci/sinpsi2);
/*  printf("f2 part3b %g \n", part3);  /* */
  part3 /= (muclabda +cosa*muci/sinpsi2)*(cosa*muclabda/sinpsi1 - muclabda);
  part3 *= cosa*exp(-(muclabda/sinpsi1 - muclabda/cosa)*thickness);
 /* printf("f2 part3c %g muclabda %g muci %g cosa %g sinpsi2 %g\n", 
	  part3, muclabda, muci,  cosa, sinpsi2);   /*  */



  result=(cross)*(part1+ part2 +part3);
/*  printf("in f2 result is %g\n", result); /* */
 /* printf(" in f2 muci is %g muclabda is %g\n", muci, muclabda); /* */ 
  f2cnt++ ;  /* registratie functiecalls */
  return(result);
}

