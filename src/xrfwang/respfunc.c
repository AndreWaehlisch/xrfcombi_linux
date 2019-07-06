/*********************************************************
*
*    COPYRIGHT (C) M.BOS 1999 , for details see copying
*
*
*   respfunc.c
*
*   contains definition for func()
*   that calculates response ( sum of squared errors)
*   for set of input parameters in structure  by
*   calculating rik's for the composition and thicknesses
*   given in the parameters and using the measured rik's
*   response is placed in structure element  par.val
*
*   M.Bos
*   Oct 1999 adjusted for Pik factors use
*
*************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

int func(struct pstruct *parin)
{
   int i, j  ;
   layer *ptop ;
   char el[3], lin[6] ;
   int z_filt;
   double mth_filt, kev;
   double totint, labda_in, primint, secint, ssq ;
   double sumfracts;  /* sum of fractions of compounds */
   extern int nr_of_meas_lines ;
   extern multicomplayer compsample ;
   extern multilayer sample; 
   extern struct spectrum_pair **pntrspectra ;
   extern meas_line lines_to_meas[MAX_LINES] ;
   extern elem_list lijst_els;
   extern double *rik_meas;
   extern double *rik_calc;
   extern double *intprep_calc;
   extern double *pi_facts;
   extern int nr_free_vars;
   extern char anode_el[4];
   extern double take_off, berryliumd ;
   extern int ref_line;

  /* test for negative params -- if so return very large ssq */
   for (i=0; i <   nr_free_vars ; i++)
   {
      if ((*parin).parm[i] < 0.0 )     /* || ((*parin).parm[i] > 1.0)) */
      {
         /* debugging 
          printf("Negative parameters found \n"); /* */
         (*parin).val=1.0e6;
         return(1);
      }
   }
 /* test for param s > 1.0 -- if so return very large ssq */

   for (i=0; i <   nr_free_vars ; i++)
   {
      if ((*parin).parm[i] > 1.0 )     /* || ((*parin).parm[i] > 1.0)) */
      {
         /* debugging  
         printf("parameters > 1.0  found \n");  /* */
         (*parin).val=1.0e5;
         return(1);
      }
   }
  
  
   
/*   printf("Na testen op neg params in func\n");  */

   put_params(&compsample, (*parin).parm);   /* put params in datastruct smpl */

/*   printf("Na put_params in respfunc \n"); /* */
   sumfracts=convcomp(&compsample, &sample);  /* to sample data (elements) for calns */
    /*  printf("sumfracts is %g \n", sumfracts); /*  */
/*   printf("Na convcomp in func \n"); /* */
   ptop=get_layer(&sample, 0);  /* for now only top layers !!! later more */
/*   printf("In func net voor ssq=0.0 \n");  /* */
   ssq=(1.0-sumfracts)*(1.0-sumfracts) ; /* later splitsen nu sum-100 */ 
/*   ssq=0.0 ; /* */
   for (j=0; j < nr_of_meas_lines ; j++)
   {
     strcpy(el, get_asymb((lines_to_meas[j]).z , &lijst_els));
     strcpy(lin, (lines_to_meas[j]).line);
     z_filt=(lines_to_meas[j]).filt_z;
     mth_filt=(lines_to_meas[j]).mth_filt; 
     kev=(lines_to_meas[j]).kev;
     totint=0.0;
/*     printf("Voor de loop over spectrum in func\n"); /* */
     for (i=0; i < PNTS_SPECTRUM ; i++ )
     { 
       labda_in=(pntrspectra[j][i]).lambda;
/*       printf("na labda_in \n"); /* */
       primint = pilabda(j, i, ptop)*(pntrspectra[j][i]).intens;
     /*  printf("na pilabda primint is %g \n",primint); /* */
       secint= intrasecflu(el, lin, labda_in, &sample, 0, primint,
                            j, i);
/*       printf("na intrasecflu  secint is %g \n", secint); /* */
       totint+=primint + secint;
/*       printf("einde lus spectrum\n"); /* */
     }
/*     printf("na afwerken alle golflengtes in func \n"); /* */
     intprep_calc[j]=totint ;
   }
   for (j=0; j < nr_of_meas_lines ; j++)
   {
       rik_calc[j]=intprep_calc[j]*pi_facts[j]/
       (intprep_calc[ref_line]*pi_facts[ref_line]);
/*     printf("rik_cal is %g \n", rik_calc[j]);  /* */
     ssq += (rik_calc[j] - rik_meas[j])*(rik_calc[j]-rik_meas[j]); /* */
/*     ssq += (1.0 - rik_calc[j]/rik_meas[j])*(1.0-
             rik_calc[j]/rik_meas[j]); */
   }
   (*parin).val=ssq ;
   /* debugging
   printf("Debugging respfunc net voor return ssq %g \n",ssq); /* */
   free_multil(&sample);
   return(0);  /* voorlopig altijd OK later ook ERROR */
}
