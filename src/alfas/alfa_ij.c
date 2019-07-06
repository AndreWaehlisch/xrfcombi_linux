/*******************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*   alfa_ij.c
*
*   calculates alfa_ij according to 
*   R.M. Rousseau and M. Bouchard
*   X-Ray Spectrom. 15 (1986) 207
*
*   coded by M.Bos
*   april 1998
*
***********************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern elem_list lijst_els;
extern meas_line lines_to_meas[MAX_LINES];
extern double psi1, psi2;
extern emudata abstabel ;

double calc_alfa_ij(multilayer *sample, int nr_line, int zi, int zj)
{
   int i, l_limit ;
   double beta_ij, mui_star, muj_star ;
   double lambda, lambdai, labd0, l_intvl ;
   double denominator, numerator;
   char linesymb[6], elsymb[3]; 

   /* calculate lambda of the measured line */
   strcpy(linesymb,"\0\0\0\0\0\0");
   strcpy(elsymb,"\0\0\0");
   strcpy(linesymb, (lines_to_meas[nr_line]).line);
   strcpy(elsymb, get_asymb((lines_to_meas[nr_line]).z, &lijst_els));
   lambdai=characwl(linesymb, elsymb, &lijst_els);
   /* debugging 
   printf("Line %s El %s lambda %g Energy %g \n",
            linesymb, elsymb, lambdai, 12.396/lambdai); /* */
   denominator=0.0;
   numerator=0.0;
   labd0=12.396/((lines_to_meas[nr_line]).kev);
   l_intvl= (LAMBDA_HIGHEST - labd0)/PNTS_SPECTRUM;
   l_limit= (lambdai-labd0)/l_intvl;
   for (i=0; i < l_limit ; i++)
   {
     /* calculate current lambda */
     lambda= labd0 + i*l_intvl;
     mui_star = calc_muilambda(zi, 12.396/lambda, &abstabel)/sin(psi1);
     mui_star += calc_muilambda(zi, 12.396/lambdai, &abstabel)/sin(psi2);
     muj_star = calc_muilambda(zj, 12.396/lambda, &abstabel)/sin(psi1);
     muj_star += calc_muilambda(zj, 12.396/lambdai, &abstabel)/sin(psi2);
     beta_ij=(muj_star/mui_star) -1.0;
     denominator += wilabdak(sample, zi, i, nr_line)*beta_ij;
     numerator += wilabdak(sample, zi, i, nr_line);
   }
   return(denominator/numerator);
}
