/****************************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*    file delta_ij.c
*
*    calculates delta_ij(labda_k) part
*    
*    cf R.M.Rousseau M. Bouchard,
*    X-Ray Spectrom. 15 (1986) 207
*
*    M.Bos
*    April 1998
******************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern elem_line *enhances[MAX_ENH_LINES];
extern meas_line lines_to_meas[MAX_LINES];
extern elem_list lijst_els;
extern emudata abstabel;

double delta_ij(int zi, int zj, int nr_line, int lambda_index, multilayer *sample)
{
   int i ;
   char el_i[3];
   double labd0, lambda, l_intvl;
   char line_i[6], line[6] ;
   elem_line *pline;
   char el_enh[3];
   layer *ptop ;
   int cnt_enhances, z_enh;
   double result, tussen, lambdaj, lambdai ;
  
   /* calculate current lambda */
   labd0=12.396/((lines_to_meas[nr_line]).kev);
   l_intvl= (LAMBDA_HIGHEST - labd0)/PNTS_SPECTRUM;
   lambda= labd0 + lambda_index*l_intvl;

   /* calculate lambda of the measured line */
   strcpy(line_i,"\0\0\0\0\0\0");
   strcpy(el_i,"\0\0\0");
   strcpy(line_i, (lines_to_meas[nr_line]).line);
   strcpy(el_i, get_asymb((lines_to_meas[nr_line]).z, &lijst_els));
   lambdai=characwl(line_i, el_i, &lijst_els);
   /* debugging  
   printf("In delta_ij Line %s El %s lambda %g Energy %g \n",
            line_i, el_i, lambdai, 12.396/lambdai); /* */

   ptop = sample->lagen[0];
   cnt_enhances = find_enhanc(el_i, line_i, ptop, lambda);
   /* debugging 
   printf("cnt_enhances in delta_ij is %d \n", cnt_enhances); /* */
   result=0.0;

   for (i=0; i < cnt_enhances; i++)
   {
     pline = enhances[i];
     z_enh = (*pline).z ;
     strcpy(line, (*pline).line);
     strcpy(el_enh, get_asymb(z_enh, &lijst_els));  
     lambdaj=characwl(line, el_enh, &lijst_els);
/*     if (( z_enh == zj)&& (lambdaj > lambda)&&( lambdaj > lambdai )) */
     if (z_enh == zj)
     {
        /* debugging  
        printf(" Enhanced by line %s of %s at lambda %g \n", line, el_enh, lambdaj); /* */
        tussen = 0.5*k_j(zj, zi, nr_line, line,  lambda_index);
        /* debugging 
        dbgres1=2.0*tussen; /* */
        /* debugging 
        dbgres2=p_ij(zj, zi, nr_line, line, lambda_index, sample); 
        dbgres4= calc_muilambda(zi, 12.396/lambdaj, &abstabel);
        dbgres5= calc_muilambda(zj, 12.396/lambda, &abstabel);
        dbgres6= calc_muilambda(zi, 12.396/lambda, &abstabel); /* */
        tussen *= p_ij(zj, zi, nr_line, line, lambda_index, sample);
        tussen *=  calc_muilambda(zi, 12.396/lambdaj, &abstabel);
        tussen *= calc_muilambda(zj, 12.396/lambda, &abstabel);
        tussen /= calc_muilambda(zi, 12.396/lambda, &abstabel);
        result += tussen;
        /* debugging 
        printf("k_j %g, p_ij %g, muij %g, mujl %g, muil\n",
               dbgres1, dbgres2, dbgres4, dbgres5, dbgres6); /* */
     }
   }
   clean_up_enh(cnt_enhances);
   return(result);
}

