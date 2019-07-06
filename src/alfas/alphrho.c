/**********************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*      file alphrho.c
*      calculate all alpha_ij's and rho_ij's
*      cf R.M. Rousseau and M. Bouchard
*      X-Ray Spectrom. 15 (1986) 207
*
*      M.Bos
*      April 1998
*      
*      N.B. spectra are supposed to be present and 
*      accessable via pntrspectra !
*
***********************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern elem_list lijst_els;
extern meas_line lines_to_meas[MAX_LINES];
extern struct spectrum_pair **pntrspectra;
extern double **alpha_ij, **rho_ij ;
extern char anode_el[3];
extern double take_off, berryliumd;

void  alphrho(multilayer *sample, int nr_lines)
{
  int i, j, k ;
  layer *ptop;
  int nr_elements, z_meas, z_filt, z_elem_in_laag ;
  char eli[3], lin[6], elj[3] ;
  double mth_filt, kev ;

  ptop= sample->lagen[0];
  nr_elements=ptop->nr_elements;
  alpha_ij=(double **) malloc(nr_lines*sizeof(double *));
  for (i=0; i < nr_lines; i++)
  {
    alpha_ij[i]= (double *)malloc(nr_elements*sizeof(double));
    if (alpha_ij[i]==NULL)
    {
       printf("Cannot malloc alpha_ij\n");
       exit(1);
    }
  }
  rho_ij = (double **) malloc(nr_lines*sizeof(double *));
  for (i=0; i < nr_lines ; i++)
  {
     rho_ij[i] = (double *) malloc(nr_elements*sizeof(double ));
     if (rho_ij[i]==NULL)
     {
       printf("Cannot malloc rho_ij\n");
       exit(1);
     }
  }
  
  /* debugging 
  printf("After mallocing alpha and rho\n"); /* */
 
  for (k=0; k < nr_lines; k++)
  {
     z_meas=(lines_to_meas[k]).z;
     strcpy(eli, get_asymb(z_meas, &lijst_els));
     strcpy(lin, (lines_to_meas[k]).line);
     z_filt = (lines_to_meas[k]).filt_z;
     mth_filt = (lines_to_meas[k]).mth_filt;
     kev = (lines_to_meas[k]).kev ;

     /* debugging 
     printf("In alphrho z_filt %d, mth_filt %g, eli %s, lin %s\n", 
             z_filt, mth_filt, eli, lin); /* */
/*      pntrspectra[k] = gen_spec(anode_el, take_off, berryliumd, kev, 
                              PNTS_SPECTRUM, z_filt, mth_filt);  */
    for (j=0; j < nr_elements; j++)
    {
      z_elem_in_laag=get_z(((*ptop).elementen[j]).name, &lijst_els);
      /* debugging 
      printf("In alphrho z_elem_in_laag is %d\n", z_elem_in_laag); /* */
      if (z_elem_in_laag != z_meas)    
      {
        strcpy(elj,  ((*ptop).elementen[j]).name);
        alpha_ij[k][j] = calc_alfa_ij(sample, k, z_meas, z_elem_in_laag);
        /* debugging 
        printf("alpha_ij[%d][%d] is %g \n", k, j, alpha_ij[k][j]); /* */
        rho_ij[k][j] = calc_rho_ij(z_meas, z_elem_in_laag, k, sample);  
        /* debugging 
        printf("rho_ij[%d][%d] is %g \n", k, j, rho_ij[k][j]); /* */      
      }
    }
  }
}
