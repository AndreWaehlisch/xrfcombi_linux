/***************************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*      file pij_labda.c
*
*      calculation of Pij(labdak)
*      cf R.M.Rousseau, M. Bouchard
*      X-Ray Spectrom. 15 (1986) 207
*
*      M.Bos
*      April 1998
*
*****************************************************************/
#include  <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern elem_list lijst_els;
extern meas_line lines_to_meas[MAX_LINES];
extern double psi1, psi2 ;
extern emudata abstabel;


double p_ij(int zj, int zi, int line_nr_i, char linesymb_j[], int lambda_index, multilayer *sample)
{
   double labd0, lambda, lambdai, lambdaj, energy_in ;
   char el_j[3], linesymb[6], elsymb[3];
   double Ei, Ej;
   double mus_prime_lk, mus_lj, mus_dprime_li;  
   layer *ptop ;
   int nr_elements  ;
   double result, l_intvl;
  
   /* calculate current lambda */
   labd0=12.396/((lines_to_meas[line_nr_i]).kev);
   l_intvl= (LAMBDA_HIGHEST - labd0)/PNTS_SPECTRUM;
   lambda= labd0 + lambda_index*l_intvl;
   energy_in= 12.396/lambda;

   /* calculate lambda of the measured line */
   strcpy(linesymb,"\0\0\0\0\0\0");
   strcpy(elsymb,"\0\0\0");
   strcpy(linesymb, (lines_to_meas[line_nr_i]).line);
   strcpy(elsymb, get_asymb((lines_to_meas[line_nr_i]).z, &lijst_els));
   lambdai=characwl(linesymb, elsymb, &lijst_els);
   Ei= 12.396/lambdai;
   /* debugging 
   printf("Line %s El %s lambda %g Energy %g \n",
            linesymb, elsymb, lambdai, Ei); /* */
   
   /* calculate lambda of the enhancing line */
   strcpy(el_j,"\0\0\0");
   strcpy(el_j, get_asymb(zj, &lijst_els));
   lambdaj=characwl(linesymb_j, el_j, &lijst_els);
   Ej= 12.396/lambdaj;

   /* debugging 
   printf("Enhancing line %s El %s lambda %g Energy %g \n",
            linesymb_j, el_j, lambdaj, Ej); /* */

   ptop = (*sample).lagen[0];
   nr_elements= (*ptop).nr_elements;
 
   /* debugging 
   printf("Before calculations of mus \n"); /* */


   mus_prime_lk = calc_mumix(energy_in, ptop, &lijst_els, &abstabel)/sin(psi1);

   /* debugging 
   printf("mus_prime_lk is %g\n", mus_prime_lk); /* */
 
   mus_dprime_li = calc_mumix(Ei, ptop, &lijst_els, &abstabel)/sin(psi2);

   /* debugging 
   printf("mus_dprime_li is %g \n", mus_dprime_li); /* */

   mus_lj = calc_mumix(Ej, ptop, &lijst_els, &abstabel);

   /* debugging 
   printf("mus_lj is %g \n", mus_lj); /* */

   result =  (1.0/mus_prime_lk)*log(1.0 + mus_prime_lk/mus_lj);
   result += (1.0/mus_dprime_li)*log(1.0 + mus_dprime_li/mus_lj);
  
   return(result);
}

