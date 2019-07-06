/*********************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*   func_el.c
*
*   contains definition for func_el()
*   that calculates response ( sum of squared errors)
*   for set of input parameters in structure  by
*   calculating rxi's for the composition and thicknesses
*   given in the parameters and using the measured rxi's
*
*   M.Bos
*   Feb 1998
*
*************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

double func_el(multilayer *sample)
{
   int i, j  ;
   layer *ptop ;
   char el[3], lin[6] ;
   int z_filt;
   double mth_filt, kev;
   double totint, labda_in, primint, secint, ssq ;
   extern int nr_of_meas_lines ;
   extern struct spectrum_pair **pntrspectra ;
   extern meas_line lines_to_meas[MAX_LINES] ;
   extern elem_list lijst_els;
   extern double *rxi_meas;
   extern double *rxi_calc;
   extern double *intprep_calc;
   extern double *intbulk_calc;
   extern int nr_free_vars;
   extern char anode_el[4];
   extern double take_off, berryliumd ;

   ptop=get_layer(sample, 0);  /* for now only top layers !!! later more */
/*   printf("In func net voor ssq=0.0 \n");  /* */
   ssq=0.0 ; /* */
   for (j=0; j < nr_of_meas_lines ; j++)
   {
     strcpy(el, get_asymb((lines_to_meas[j]).z , &lijst_els));
     strcpy(lin, (lines_to_meas[j]).line);
     z_filt=(lines_to_meas[j]).filt_z;
     mth_filt=(lines_to_meas[j]).mth_filt; 
     kev=(lines_to_meas[j]).kev;
     totint=0.0;
/*     pntrspectra[j]=gen_spec(anode_el, take_off, berryliumd, kev,
                            PNTS_SPECTRUM, z_filt, mth_filt); */
/*     printf("Voor de loop over spectrum in func\n"); /* */
     for (i=0; i < PNTS_SPECTRUM ; i++ )
     { 
       labda_in=(pntrspectra[j][i]).lambda;
/*       printf("na labda_in \n"); /* */
       primint = pilabda(j, i, ptop)*(pntrspectra[j][i]).intens;
/*       printf("na pilabda primint is %g \n",primint); /* */
       secint= intrasecflu(el, lin, labda_in, sample, 0, primint,
                            j, i);
/*       printf("na intrsecflu  secint is %g \n", secint); /* */
       totint+=primint + secint;
/*       printf("einde lus spectrum\n"); /* */
     }
/*     printf("na afwerken alle golflengtes in func \n"); /* */
     intprep_calc[j]=totint ;
     rxi_calc[j]=totint/intbulk_calc[j];
     /* debugging
     printf("rxi_cal is %g meas %g\n", rxi_calc[j], rxi_meas[j]);  /* */
     ssq += (rxi_calc[j] - rxi_meas[j])*(rxi_calc[j]-rxi_meas[j]); /* */
  ssq += (1.0 - rxi_calc[j]/rxi_meas[j])*(1.0 - rxi_calc[j]/rxi_meas[j]); /* */
   }
   /* debugging  
   printf("Debugging respfunc net voor return ssq %g \n",ssq); /* */
   return(ssq);  
}
