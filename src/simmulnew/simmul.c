/***********************************************************************
*
*                    COPYRIGHT (C) M.Bos, FOR DETAILS SEE FILE COPYING
*
*                    simmul.c
*                    Program for quantitative 
*                    X-ray fluorescence analysis
*                    in multilayer samples 
*                    to calculate example from
*                    D.K.G. de Boer, X-Ray spectrom. 19 (1990) 145
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                    september 1998
*
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

int  main(int argc, char *argv[])
{
   multilayer sample ;
   multilayer bulksample ;
   layer  *pbulk, *pcurrlayer, *potherlayer ;
   double primint, secint, totint, totint_bulk;
   struct spectrum_pair *pntrspectrum ;
   double labda_in, labda_i, energy_in, energy_i ;
   int i, j, k, l, nr_layers ;
   char el[3];
   char lin[6];
   int z_filt;
   double mth_filt, kev;
   FILE *fpresult;
   double currthickn, otherthickn ;
   double mualprim, muaiprim ; 
   double corr_ins ;
    
   fpresult=fopen("uitvoer","w");
   init_tables();
   /* debugging
   printf("na init tables \n"); /* */

   set_spectrom();
   /* debugging
   printf("na set spectrom \n"); 
   for (i=0; i < nr_of_meas_lines; i++)
   {
       printf("%d Z: %d lijn %s HT %f mthickn filt %g  \n",
               i, (lines_to_meas[i]).z, (lines_to_meas[i]).line,
                (lines_to_meas[i]).kev, (lines_to_meas[i]).mth_filt);
   } /* */
   def_sample(&sample);
   checkdata(&sample);  /* debugging pa0mbo */
   nr_layers = sample.nr_layers ;


   for (j=0; j < nr_of_meas_lines; j++) /* <========================= */
   {
     strcpy(el, get_asymb((lines_to_meas[j]).z, &lijst_els));
     strcpy(lin, (lines_to_meas[j]).line);
 
     /* debugging 
     printf("in main lijn %s element Z %d %s el %s\n", lin,
           (lines_to_meas[j]).z,
             get_asymb((lines_to_meas[j]).z, &lijst_els),el); /* */

     z_filt=(lines_to_meas[j]).filt_z;
     mth_filt=(lines_to_meas[j]).mth_filt;

     /* debugging 
     printf("Filter Z %d, filter massthickness %g \n",
            z_filt, mth_filt); /* */
 
     kev=(lines_to_meas[j]).kev;
     pntrspectrum=gen_spec(anode_el, take_off, berryliumd, kev,
                 PNTS_SPECTRUM, z_filt, mth_filt);
     labda_i = characwl(lin, el, &lijst_els);
     energy_i= 12.396/labda_i ;
     pbulk=mkbulklay(el);
     /* now set up complete sample for bulk element */
      bulksample.nr_layers= 1;
      bulksample.lagen=&pbulk; 
     totint=0.0;
     totint_bulk=0.0 ;
      
      for (i=0 ; i < PNTS_SPECTRUM ; i++)  /* <============== */
      {
        labda_in= (pntrspectrum[i]).lambda ;
        energy_in = 12.396 / labda_in ;
/* debugging
 printf("Voor pilamba i is %d en intens %g \n",i, (pntrspectrum[i]).intens); /* */ 
      totint_bulk += pilabda(el,lin,labda_in, pbulk)*(pntrspectrum[i]).intens;
      for (k=0; k < nr_layers ; k++)   /* <================================ */
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
      /* debugging 
   printf("totint %s %s is %e bulkint %e\n",el, lin, totint,totint_bulk); /*  */
     printf("Rel. int %s %s is %e \n",el,lin, totint/totint_bulk);
     fprintf(fpresult,"%s %s %e\n", el, lin, totint/totint_bulk);
 }
  fclose(fpresult);
  return 0;
}

