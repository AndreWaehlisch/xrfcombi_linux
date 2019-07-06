/***********************************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*                    xrf.c
*                    Main program for quantitative 
*                    X-ray fluorescence analysis
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                    November 1997
*
*		     update 8 jan 98
*                    store generated spectra for later use
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

#define MAXITERS 800
double *free_pars ;  /* pointer to array with simplex params + resp */
double *pars_min ;
double *pars_max; 

int main(int argc,char *argv[])
{
   extern multicomplayer compsample ;   /* sample in terms of compounds */
   extern multilayer sample ;          /* sample in terms of elements */
   extern struct spectrum_pair **pntrspectra ;
   FILE *fp;
   char filename[80]="sample.dat"; 
   char dirname[80]="./data/";
   extern char title[132]; 
   extern double *rxi_meas; 
   extern double *rxi_calc;
   extern double *intbulk_calc;
   extern double *intprep_calc;
   extern struct pstruct *p, pcent, **p_p ;
   extern int prt_cycle;
   extern int ndata;
   extern double quad_test, test, ypmin, yzero; 
   extern int quad_cycle, nquad_skip, maxquad_skip;
   double thickness, old_thickn;
   int itercnt ;
/* debugging
   printf("voor init tables \n"); /* */
   init_tables();
/* debugging
   printf("na init tables \n");  /* */

   set_spectrom();
   /* debugging
   printf("na set spectrom \n");
   printf("nr_of_meas_lines %d \n",nr_of_meas_lines); /* */

   rxi_meas= malloc(nr_of_meas_lines* sizeof(double));
   rxi_calc= malloc(nr_of_meas_lines* sizeof(double));
   intbulk_calc= malloc(nr_of_meas_lines* sizeof(double));
   intprep_calc= malloc(nr_of_meas_lines* sizeof(double));
   
   fp= own_fopen("./data/lines.dat");
   get_rximeas(fp, nr_of_meas_lines, rxi_meas);
   fclose(fp);
   /* debugging
   printf("nal lezen rxi's \n"); 
   for (i=0; i < nr_of_meas_lines; i++)
   {
       printf("%d Z: %d lijn %s rel int. %g \n", i, (lines_to_meas[i]).z,
              (lines_to_meas[i]).line, rxi_meas[i]);
   }  /* */


   /* allocate memory for spectra per measured line */
   pntrspectra= (struct spectrum_pair **) malloc( nr_of_meas_lines*
                 sizeof(struct spectrum_pair *));


/*   printf("Give filename with sample data ");
   gets(filename) ;  */
   strcpy(dirname+7,filename);
   fp=own_fopen(dirname);
   bldcompos(fp, &compsample);
   fclose(fp);
   convcomp(&compsample, &sample);

   /* now set up and fill block with abscoeff for elemens with labda */
   /* debugging
   printf("Voor call fillmulabda\n"); /* */
   fillmulabda(&sample,  nr_of_meas_lines);
   /* debugging
   printf("Na call fillmulabda\n");  /* */

   /* Now only once for bulk materials */

   calc_alphas(&sample, nr_of_meas_lines);
   /* debug 
   printf("Na calc_alphas\n");  /* */
   conciter(&sample); 
/* debugging  
   printf("na conciter \n"); /* */
   old_thickn=1.0e3;
   thickness=((&sample)->lagen[0])->massthickness;
   /* debugging 
   printf("thickness layer 0 is %g \n",thickness); /* */
   itercnt=100; 
   while ((fabs(1.0-thickness/old_thickn) > 1.0e-4)&& (itercnt >0))
   {
     criss(&sample);
     /* debugging  */
     printf("sum of fracts el is %g \n",sum_el(&sample,0)); /* */
     old_thickn=thickness;
     thickness= thickness*sum_el(&sample,0);
     /* debugging */
     printf("curr. thickness is %g\n",thickness);  /* */
     ((&sample)->lagen[0])->massthickness=thickness;
     itercnt--;
   }
   return(0);
}
