/*********************************************************
*
*    COPYRIGHT (C) M.BOS 2004 , for details see copying
*
*
*   respfuncscat.c
*
*   contains definition for funcscat()
*   that calculates response ( sum of squared errors)
*   for set of input parameters in structure  by
*   calculating rxi's for the composition and thicknesses
*   given in the parameters and using the measured rxi's
*   response is placed in structure element  par.val
*
*   M.Bos
*   April 12th 2004
*
*************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

int funcscat(struct pstruct *parin)
{
   int i, j  ;
   multilayer bulksample;
   layer *ptop, *pbulk ;
   char el[3], lin[6] ;
   int z_filt;
   double mth_filt, kev;
   double  labda_in,  ssq ;
   double primint, secint, totint, comptonint, primint_bulk, comptonint_bulk, totint_bulk;
   double  totcomptonint, totcomptonint_bulk, totprimint, totprimint_bulk ;
   double totsecint ;
   double raylfotint, raylfotint_bulk, totraylfotint, totraylfotint_bulk;
   double fotraylint, fotraylint_bulk, totfotraylint, totfotraylint_bulk;
   double sumfracts;  /* sum of fractions of compounds */
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

  /* test for negative params -- if so return very large ssq */
   for (i=0; i <   nr_free_vars ; i++)
   {
      if ((*parin).parm[i] < 0.0 )     /* || ((*parin).parm[i] > 1.0)) */
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
      if ((*parin).parm[i] > 1.0 )     /* || ((*parin).parm[i] > 1.0)) */
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
   sumfracts=convcomp(&compsample, &sample);  /* to sample data (elements) for calns */
/*     printf("sumfracts is %g \n", sumfracts); */
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
     pbulk=mkbulklay(el);
     bulksample.nr_layers=1;
     bulksample.lagen=&pbulk;

     totint=0.0;
     totsecint=0.0;
     totint_bulk=0.0;
     totcomptonint=0.0;
     totcomptonint_bulk=0.0;
     comptonint=0.0 ;
     comptonint_bulk=0.0 ;
     totprimint=0.0;
     totprimint_bulk=0.0;
     totraylfotint =0.0;
     totraylfotint_bulk=0.0;
     raylfotint=0.0;
     raylfotint_bulk =0.0;
     fotraylint =0.0;
     fotraylint_bulk = 0.0;
     totfotraylint = 0.0;
     totfotraylint_bulk = 0.0;
   
     for (i=0; i < PNTS_SPECTRUM ; i++ )
     { 
       labda_in=(pntrspectra[j][i]).lambda;
      primint= pilabdascat(el,lin,labda_in, ptop)*(pntrspectra[j][i]).intens;
      primint_bulk = pilabdascat(el,lin,labda_in, pbulk)*(pntrspectra[j][i]).intens;
      secint=intrasecflu_scat(el,lin,labda_in, &sample, 0, primint) ; 
      comptonint=comptisesr(el,lin, labda_in, ptop,i, 0) *(pntrspectra[j][i]).intens; 
      comptonint_bulk=comptisesr(el, lin,labda_in, pbulk,i, 1)*(pntrspectra[j][i]).intens; 
      raylfotint=raylisesr(el,lin, labda_in, ptop,i, 0) *(pntrspectra[j][i]).intens; 
      raylfotint_bulk=raylisesr(el, lin,labda_in, pbulk,i, 1)*(pntrspectra[j][i]).intens; 
      fotraylint=ispfr(el,lin, labda_in, ptop,i, 0) *(pntrspectra[j][i]).intens; 
      fotraylint_bulk=ispfr(el, lin,labda_in, pbulk,i, 1)*(pntrspectra[j][i]).intens; 
      /* debugging  
      printf("pnt %d %g %g %g %g %g %g %g %g %g\n", i,primint, primint_bulk,
		  secint, comptonint, comptonint_bulk, raylfotint, raylfotint_bulk,
		  fotraylint, fotraylint_bulk); /* */

      totprimint += primint;
      totsecint += secint;
      totcomptonint +=comptonint;
      totraylfotint += raylfotint;
      totfotraylint += fotraylint;
      totint += primint+secint+comptonint+raylfotint+ fotraylint;

      totprimint_bulk += primint_bulk; 
      totcomptonint_bulk+=comptonint_bulk;
      totraylfotint_bulk+=raylfotint_bulk;
      totfotraylint_bulk+=fotraylint_bulk;
      totint_bulk +=primint_bulk+comptonint_bulk+raylfotint_bulk+fotraylint_bulk;
     }
  /*  printf("%s %s totint %g totint_bulk %g totprimint %g totprimintbulk %g"
                     "totsecint %g totcomptonint %g totcomptonintbulk %g\n"
               "totraylfotint %g totraylfotint_bulk %g \n"
               "totfotraylint %g totfotraylint_bulk %g  \n",
            el, lin, totint, totint_bulk, totprimint, totprimint_bulk,
                         totsecint, totcomptonint, totcomptonint_bulk,
                  totraylfotint, totraylfotint_bulk, totfotraylint, totfotraylint_bulk);
	printf("rxi old is %g rxi scatter is %g \n", 
		(totprimint+totsecint)/(totint_bulk-totcomptonint_bulk-totraylfotint_bulk-
		   totfotraylint_bulk), totint/totint_bulk); /* */
    rxi_calc[j]=totint/totint_bulk;

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
