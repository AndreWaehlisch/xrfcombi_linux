/*********************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*   respfunc.c
*
*   contains definition for func()
*   that calculates response ( sum of squared errors)
*   for set of input parameters in structure  by
*   calculating rxi's for the composition and thicknesses
*   given in the parameters and using the measured rxi's
*   response is placed in structure element  par.val
*   Now all layers are taken into account 
*
*   M.Bos
*   Sept  1998
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
   double labda_i, energy_i, energy_in ;
   int nr_layers ;
   extern int nr_of_meas_lines ;
   extern multicomplayer compsample ;
   extern multilayer sample; 
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
   extern double sumfracts[MAX_LAYERS] ;

  /* test for negative params -- if so return very large ssq */
   for (i=0; i <   nr_free_vars ; i++)
   {
      if ((*parin).parm[i] < 0.0 )  
      {
         /* debugging
          printf("Negative parameters found \n"); /* */
         (*parin).val=1.0e4;
         return(1);
      }
   }
 /* test for param s > 1.0 -- if so return very large ssq */

   for (i=0; i <   nr_free_vars ; i++)
   {
      if ((*parin).parm[i] > 1.0 )     
      {
         /* debugging
         printf("parameters > 1.0  found \n");  /* */
         (*parin).val=1.0e4;
         return(1);
      }
   }
  
  
   
/*   printf("Na testen op neg params in func\n");  */

   put_params(&compsample, (*parin).parm);   /* put params in datastruct smpl */

/*   printf("Na put_params in respfunc \n"); /* */
   convcomp(&compsample, &sample);  /* to sample data (elements) for calns */
/*   printf("Na convcomp in func \n"); /* */

   nr_layers = sample.nr_layers ;
   
   ssq=0.0 ; 
   for (i=0 ; i < nr_layers ; i++)
   {
    ssq += (1.0-sumfracts[i])*(1.0-sumfracts[i]) ; 
   }

   for (j=0; j < nr_of_meas_lines ; j++)
   {
     strcpy(el, get_asymb((lines_to_meas[j]).z , &lijst_els));
     strcpy(lin, (lines_to_meas[j]).line);
     z_filt=(lines_to_meas[j]).filt_z;
     mth_filt=(lines_to_meas[j]).mth_filt; 
     kev=(lines_to_meas[j]).kev;
     labda_i = characwl(lin, el, &lijst_els);
     energy_i = 12.396/labda_i ; 

     totint=0.0;

/*     printf("Voor de loop over spectrum in func\n"); /* */
     for (i=0; i < PNTS_SPECTRUM ; i++ )
     { 
       labda_in=(pntrspectra[j][i]).lambda;
/*       printf("na labda_in \n"); /* */

       energy_in = 12.396 / labda_in ;
   
       for (k=0; k < nr_layers; k++)
       { 
         pcurrlayer = get_layer(&sample, k);
         currthickn = pcurrlayer->massthickness ; 
         corr_ins = 0.0 ;
         l =0; 
         while ( l < k)
         {
           potherlayer = get_layer(&sample, l);
           otherthickn= potherlayer->massthickness ;

           /* debug 
           printf("thickness other layer is %e \n", otherthickn); /* */

           mualprim= calc_mumix(energy_in, potherlayer, &lijst_els, &abstabel);
           muaiprim= calc_mumix(energy_i, potherlayer, &lijst_els, &abstabel);

           /* debug 
           printf("mualprim is %e muaiprim is %e \n", mualprim, muaiprim);/* */

           corr_ins -= (mualprim/sin(psi1) + muaiprim/sin(psi2))*otherthickn ;
           l++ ;
         }
         corr_ins = exp(corr_ins);

         /* debugging 
         printf(" corr_ins factor is %e \n", corr_ins); /* */

         primint = pilabda(el,lin,labda_in, pcurrlayer)*
                   corr_ins*(pntrspectrum[i]).intens; 

         secint = intrasecflu(el,lin, labda_in , &sample, k, primint);

         /* debugging 
         printf("primint %g  secint %g\n",
                primint, secint);  /* */

         totint += primint+secint;

        /*debugging 
  printf("totint na intra %g totint_bulk %g \n", totint, totint_bulk); /* */

        l=0;
        while (l < k )
        {
           totint += intersecflu(el, lin, labda_in, &sample, k, l, primint);
           l++ ;
        }
        l= k+1 ;
        while ( l < nr_layers)
        {
           totint += intersecflu(el, lin, labda_in, &sample, k, l, primint); 
           l++ ;
        }

     }
   }
   intprep_calc[j] = totint/intbulk_calc[j] ;
      /* debugging */
   printf("totint %s %s is %e bulkint %e\n",el, lin, totint,intbulk_calc[j]); /* */
     printf("Rel. int %s %s is %e \n",el,lin, intprep_calc[j]);
/*     printf("na afwerken alle golflengtes in func \n"); /* */
     intprep_calc[j]=totint ;
     rxi_calc[j]=totint/intbulk_calc[j];
/*     printf("rxi_cal is %g \n", rxi_calc[j]);  /* */
     ssq += (rxi_calc[j] - rxi_meas[j])*(rxi_calc[j]-rxi_meas[j]); /* */
  ssq += (1.0 - rxi_calc[j]/rxi_meas[j])*(1.0 - rxi_calc[j]/rxi_meas[j]); /* */
   }
   (*parin).val=ssq ;
   /* debugging
   printf("Debugging respfunc net voor return ssq %g \n",ssq); /* */
   free_multil(&sample);
   return(0);  /* voorlopig altijd OK later ook ERROR */
}
