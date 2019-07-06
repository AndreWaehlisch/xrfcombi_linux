/*******************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*    wilabdak.c
*    calculates function Wi(labda_k)
*    from paper R.M.Rousseau and M.Bouchard
*    X-Ray Spectrom. 15 (1986)207
*
*********************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern elem_list lijst_els;
extern meas_line lines_to_meas[MAX_LINES];
extern double muzlambda[MAX_LINES][TOTAL][PNTS_SPECTRUM];
extern struct spectrum_pair **pntrspectra;
extern double psi1, psi2 ;
extern emudata abstabel ;

double wilabdak(multilayer *sample, int zi, int lambda_index,
                int nr_line)
{
   int i,  zj ;
   double beta_ij, mui_star, muj_star, sumcbeta ;
   int nr_elements ;
   layer *firstlayer ;
   elem *pelementen ;
   char linesymb[6];
   char elsymb[3];
   double lambda, lambdai, labd0, l_intvl;
   double result;

   /* calculate current lambda */
   labd0=12.396/((lines_to_meas[nr_line]).kev); /* belongs to kev line i */
   l_intvl= (LAMBDA_HIGHEST - labd0)/PNTS_SPECTRUM;
   lambda= labd0 + lambda_index*l_intvl;

   /* calculate lambda of the measured line */
   strcpy(linesymb,"\0\0\0\0\0\0");
   strcpy(elsymb,"\0\0\0");
   strcpy(linesymb, (lines_to_meas[nr_line]).line);
   strcpy(elsymb, get_asymb((lines_to_meas[nr_line]).z, &lijst_els));
   lambdai=characwl(linesymb, elsymb, &lijst_els);
   /* debugging 
   printf("Line %s El %s lambda %g Energy %g \n",
            linesymb, elsymb, lambdai, 12.396/lambdai); /* */
   
   firstlayer=sample->lagen[0];
   nr_elements=firstlayer->nr_elements;
   pelementen=firstlayer->elementen ;
   /* debugging 
   printf("Nr of elements in first layer is %d\n",nr_elements); /* */
   sumcbeta=0.0;
   for (i=0; i < nr_elements; i++)
   {
     zj=get_z(pelementen[i].name,lijst_els);
     /* debugging 
     printf("found element %s with Z=%d\n", pelementen[i].name,zj); /* */
     mui_star = calc_muilambda(zi, 12.396/lambda, &abstabel)/sin(psi1);
     mui_star += calc_muilambda(zi, 12.396/lambdai, &abstabel)/sin(psi2);
     if (zj !=zi)
     {
       muj_star = calc_muilambda(zj, 12.396/lambda, &abstabel)/sin(psi1);
       muj_star += calc_muilambda(zj, 12.396/lambdai, &abstabel)/sin(psi2); 
       beta_ij= muj_star/mui_star - 1.0 ;
       sumcbeta += pelementen[i].wfract*beta_ij;
     }
   }
   result = calc_muilambda(zi, 12.396/lambda, &abstabel);
   result *= pntrspectra[nr_line][lambda_index].intens*l_intvl;
   result /= mui_star;
   result /= (1.0 + sumcbeta);
   return(result);
}
   
