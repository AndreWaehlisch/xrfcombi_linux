/*
*    file comptoninteg.c
*
*    functions to calculate integrand
*    for double integral for crosssection
*    compton scattering
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
*    double store_labda
*    double store_labdai
*
*    and two dependent variables alfa and fi
*   
*    now based on paper G.Tirao and G. Stutz
*     XRS 32(2003) 13
*
*    code by M.Bos
*
*    Jan 17th 2004
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
extern double store_labdai  ;
extern emudata abstabel ;
extern int store_z ;
extern layer *ptrlaag ;
extern double store_energy;
extern double store_thickness;
extern double fluor_grens;

double comptfot_func1(double cosa, double fi)
{
  double theta, costheta, cross, result, thetadegr;
  double tauilambdastar, lambdastar , energy, energy_star ;
  double mulambdastar, alfa, part1, part2, part3 ;
  double thickness ;
  alfa = acos(cosa);

 /* printf("in func1 store_energy is %g \n", store_energy); /*  */
  /* calc theta for case t1 < t2 */
 /* printf("theta1 is %g \n",theta); /* */
  costheta= sinpsi1*cosa+sin(alfa)*cos(fi)*cospsi1;
  thickness = store_thickness;
 /* printf("costheta is %g \n",costheta); /* */
  theta= acos(costheta);
  lambdastar= store_labda + Compton*(1.0-costheta);

  /* added condition Jan 19th 2004 */
  if (lambdastar > fluor_grens)  return(0.0); /* not fluoresced */
  energy_star= 12.396/lambdastar ;
  mulambdastar = calc_mumix(energy_star, ptrlaag, &lijst_els, &abstabel);
  tauilambdastar = calc_muilambda(store_z, energy_star, &abstabel);
  
 /* printf("theta is %f \n",theta); /* */
  thetadegr= (theta/PI)*180.0;
  /* debugging */
 /* printf("in f1 theta is %g \n",theta); /* */
  /* debugging 
   printf("in f1 alfa %g fi %g costhe %g theta %g sinpsi2 %g cospsi2 %g PI %g\n", 
         alfa,fi, costheta,theta, sinpsi2, cospsi2, PI); /* */
    cross=comptcrossmix(store_energy, ptrlaag, &lijst_els, thetadegr); /* */
 /*  printf("in f1 compton cross is %g bij thetadegr %g\n", cross, thetadegr); /*  */
    part1 = 1.0/(mulambdastar+cosa*muci/sinpsi2);
/*	printf("f1 part1 %g mulambdastar %g muci %g\n", 
		     part1, mulambdastar, muci); /* */
	part2 = exp(-(muclabda/sinpsi1 + muci/sinpsi2)*thickness);
	part2 /= (cosa*muclabda/sinpsi1 - mulambdastar);
/*	printf("f1 part2 %g muclabda %g mulambdastar %g cosa  %g \n", 
		    part2, muclabda, mulambdastar, cosa); /*  */
	part3 = -cosa*exp(-(mulambdastar/cosa+muci/sinpsi2)*thickness);
	part3 *= (muclabda/sinpsi1+muci/sinpsi2);
	part3 /= (mulambdastar + cosa*muci/sinpsi2)*(cosa*muclabda/sinpsi1-mulambdastar);
/*    printf("f1 part3 %g mulambdastar %g muci %g cosa %g \n", 
		   part3, mulambdastar, muci, cosa);  /* */

    result=(cross*tauilambdastar)*(part1+part2+part3); 
/*	printf("in f1 result is %g at alfa %g and fi %g thetadegr %g\n", result, alfa, fi, thetadegr);  /* */
  return(result);
}


double comptfot_func2(double cosa, double fi)
{
  double theta, costheta, cross, result, thetadegr;
  double tauilambda, lambdastar , energy, energy_star ;
  double mulambdastar, lambda, alfa,  thickness, part1, part2, part3 ;
  double tauilambdastar ;
  /* calc theta for case t1 > t2 */
 
   alfa=acos(cosa);

  costheta= sinpsi1*cos(alfa)-sin(alfa)*cos(fi)*cospsi1;
  theta= acos(costheta);
  thetadegr= (theta/PI)*180.0 ;
  thickness = store_thickness;
  lambdastar= store_labda + Compton*(1.0-costheta);
  /* added condition Jan 19th 2004 */
  if (lambdastar > fluor_grens)  return(0.0); /* not fluoresced */
  energy= 12.396/store_labda ;
  energy_star = 12.396/lambdastar ;
  mulambdastar = calc_mumix(energy_star, ptrlaag, &lijst_els, &abstabel);
  tauilambdastar = calc_muilambda(store_z, energy_star, &abstabel);
  /* debugging 
  printf("in f2 alfa %g fi %g theta %g sinpsi2 %g cospsi1 %g \n", 
         alfa,fi,theta, sinpsi2, cospsi1); /* */
  /* cross=0.1234 ; /* ter test cross opzoeken uitgeschakeld */
   cross=comptcrossmix(store_energy, ptrlaag, &lijst_els, thetadegr); /*  */
 /*  printf("in f2  alfa is %g  cos %g tan %g , fi %g \n",alfa,  cos(alfa), tan(alfa), fi); /* */
  /* printf("in f2 cross is %g bij theta %g\n,", cross, thetadegr); /* */
    part1 = 1.0/(mulambdastar-cosa*muclabda/sinpsi1);
/*	printf("f2 part1 %g mulambdastar %g cosa %g muclabda %g \n",
		     part1, mulambdastar, cosa, muclabda);  /*  */
	part2 = exp(-(muclabda/sinpsi1+muci/sinpsi2)*thickness);
	part2 /= (-cosa*muci/sinpsi2 - mulambdastar);
/*	printf("f2 part2 %g muclabda %g mulambdastar %g cosa %g\n",
		   part2, muclabda, mulambdastar, cosa);  /*  */
	if (cosa==0)
		part3 = 0.0;
	else
    {

	    part3= cosa*exp(-(muclabda/sinpsi1-mulambdastar/cosa)*thickness);
	    part3*= (muclabda/sinpsi1 + muci/sinpsi2);
	    part3 /= (mulambdastar+cosa*muci/sinpsi2)*(cosa*muclabda/sinpsi1-mulambdastar);
	}
  /*  printf("f2 part3 %g cosa %g muclabda %g muci %g mulambdastar %g\n",
		    part3, cosa,  muclabda, muci, mulambdastar) ;  /*  */
    

  result=cross*tauilambdastar*(part1 + part2 + part3);

/*  printf("in f2 result is %g\n", result); /* */
 /* printf(" in f2 muci is %g muclabda is %g\n", muci, muclabda); /* */ 
  f2cnt++ ;  /* registratie functiecalls */
  return(result);
}

