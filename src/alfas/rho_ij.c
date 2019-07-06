/*********************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*  file rho_ij.c
*
*  routine to calculate rho_ij
*  from paper by R.M. Rousseau and M.Bouchard
*  X-Ray Spectrom. 15 (1986) 207
*
*  M.Bos
*  April 1998
*
*********************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern meas_line lines_to_meas[MAX_LINES];
extern struct spectrum_pair **pntrspectra ; 
extern elem_list lijst_els;

double calc_rho_ij(int zi, int zj, int nr_line, multilayer *sample)
{
  int i, l_limit;
  double  denominator, numerator , lambdai, labd0, l_intvl;
  char linesymb[6], elsymb[3]; 
 
  denominator = 0.0;
  numerator = 0.0;

  /* calculate lambda of the measured line */
  strcpy(linesymb,"\0\0\0\0\0\0");
  strcpy(elsymb,"\0\0\0");
  strcpy(linesymb, (lines_to_meas[nr_line]).line);
  strcpy(elsymb, get_asymb((lines_to_meas[nr_line]).z, &lijst_els));
  lambdai=characwl(linesymb, elsymb, &lijst_els);
  /* debugging 
  printf("Line %s El %s lambda %g Energy %g \n",
           linesymb, elsymb, lambdai, 12.396/lambdai); /* */
  labd0=12.396/((lines_to_meas[nr_line]).kev);
  l_intvl= (LAMBDA_HIGHEST - labd0)/PNTS_SPECTRUM;
  l_limit= (lambdai-labd0)/l_intvl;
  for (i=0; i < l_limit ; i++)
  {
    denominator += wilabdak(sample, zi, i, nr_line)*delta_ij(zi, zj,
                   nr_line, i, sample);

    /* debugging 
    printf("denominator in calc_rho_ij is %g \n", denominator); /* */

    numerator += wilabdak(sample,  zi,  i, nr_line);
  }
  return(denominator/numerator);
}  
  
    

