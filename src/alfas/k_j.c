/********************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*    file k_j.c
*
*    calculation of k_j factor in paper by
*    R.M. Rousseau and M. Bouchard
*    X-Ray Spectrom. 15 (1986) 207
*
*    M.Bos
*    April 1998
*
***********************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern double omegas[TOTAL][NR_OMEGAS] ;
extern double relrates[TOTAL][TOT_RELRATES] ;
extern double jmpdata[ALL_ELMS][2];
extern meas_line lines_to_meas[MAX_LINES];
extern char *shellreeks[27];
extern char *linenames[27];
extern elem_list lijst_els ;
double labd0;

double k_j(int zj, int zi, int nr_line_i, char linesymb_j[], int lambda_index)
{
   double lambda, lambdai, energy_in, result, l_intvl ;
   char linesymb[6], elsymb[3];
   int hole_nr;
   char *shellname_j ;
 

   /* debugging 
   printf("On entry k_j sinemsymb_j is %s \n", linesymb_j); /* */   
   /* calculate current lambda */
   labd0=12.396/((lines_to_meas[nr_line_i]).kev);
   l_intvl= (LAMBDA_HIGHEST - labd0)/PNTS_SPECTRUM;
   lambda= labd0 + lambda_index*l_intvl;

   /* calculate lambda of the measured line */
   strcpy(linesymb,"\0\0\0\0\0\0");
   strcpy(elsymb,"\0\0\0");
   strcpy(linesymb, (lines_to_meas[nr_line_i]).line);
   strcpy(elsymb, get_asymb((lines_to_meas[nr_line_i]).z, &lijst_els));
   lambdai=characwl(linesymb, elsymb, &lijst_els);
   /* debugging 
   printf("In k_j linesymb %s elsymb %s lambdai %g \n", linesymb, elsymb, lambdai); /* */
   /* debugging 
   printf("Line %s El %s lambda %g Energy %g \n",
            linesymb, elsymb, lambdai, 12.396/lambdai); /* */
 
   energy_in=12.396/lambda;
   shellname_j  = line_to_shell(linesymb_j, linenames, shellreeks);
   hole_nr = get_hole_nr(shellname_j);
   result = getjmpfact(zj, hole_nr, energy_in, jmpdata)*getrelrate(zj, linesymb_j, relrates);
   result *= get_omega(zj, shellname_j, omegas); 
   return(result);
}


 
     

    

