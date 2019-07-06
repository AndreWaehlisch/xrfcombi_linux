/***********************************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*                    mulabfil.c
*                    routine to calculate all pure element
*                    abs coeff. at all wavelengths
*                    needed later
*                    M.BOS
*                    JAN 1998
*
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

extern elem_list lijst_els ;
extern double muzlambda[MAX_LINES][TOTAL][PNTS_SPECTRUM] ;
extern meas_line  lines_to_meas[MAX_LINES] ;

void fillmulabda(multilayer *sample, int nr_lines)
{
   int  nr_layers, tot_elements;
   layer **plagen, *plaag;
   elem *pelem_in_laag;
   int z_elem_in_laag;
   double lambda, labd0, l_intvl ;
   int i,j,k, l ;
   extern emudata abstabel;

   tot_elements=0; 
   nr_layers=(*sample).nr_layers;
   plagen= (*sample).lagen;
   for (i=0; i < nr_layers; i++)
   {
    plaag=plagen[i];
    tot_elements += (*plaag).nr_elements;
   }
   /* debugging
   printf("Totaal aantal elementen in preparaat %d \n",tot_elements); /* */
   for (k=0; k < nr_lines; k++)
   {
    labd0=12.396/((lines_to_meas[k]).kev);
    l_intvl= (LAMBDA_HIGHEST -labd0)/PNTS_SPECTRUM;
    for (i=0; i < nr_layers; i++)
    {
     plaag=plagen[i];
     pelem_in_laag= (*plaag).elementen;
     for (j=0; j< (*plaag).nr_elements ; j ++)
     {
       z_elem_in_laag= get_z(((*plaag).elementen[j]).name, &lijst_els);
       /* printf("Z-elem in muzlambda %d van el %s \n", z_elem_in_laag,
               ((*plaag).elementen[j]).name); /* */
       for (l=0; l < PNTS_SPECTRUM; l++)
       {
        lambda=labd0+l*l_intvl;
     /*    printf("In spectrumlus lambda is %g\n", lambda); /* */
        muzlambda[k][z_elem_in_laag][l]=
        calc_muilambda(z_elem_in_laag, 12.396/lambda,
                        &abstabel);
/*        printf("muzlambda[%d][%d][%d] = %g \n",
                k, z_elem_in_laag,l, muzlambda[k][z_elem_in_laag][l]); /* */

       }
     }
   }  
  }
}


void check_mu_el(int z, int line_nr)
{
   int i,j ;

   for (i=0; i < PNTS_SPECTRUM/10; i++)
   {
     for (j=0; j < 10 ; j++)
     {
      printf("%g ", muzlambda[line_nr][z][i*10 +j]);
     }
   }
}

   
