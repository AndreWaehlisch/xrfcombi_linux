/***********************************************************************
*
*    COPYRIGHT (C) M.BOS 1999 , for details see copying
*
*
*                    xrfrous2.c
*                    Main program for quantitative 
*                    X-ray fluorescence analysis
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                    September 1999
*                    Incorporated simplex to calculate
*                    compound composition after 1st
*                    round of Rousseau algorithm
*                    to obtain better estimate
*                    of not-measured elements
*
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

#define MAXITERS 1000
double *free_pars ;  /* pointer to array with simplex params + resp */
double *pars_min ;
double *pars_max; 

int main(int argc,char *argv[])
{
   extern multicomplayer compsample ;   /* sample in terms of compounds */
   extern multilayer sample ;          /* sample in terms of elements */
   extern struct spectrum_pair **pntrspectra ;
   char filename[80]="sample.dat"; 
   char dirname[80]="./data/";
   extern char title[132]; 
   extern double *rxi_meas; 
   extern double *rxi_calc;
   extern double *intbulk_calc;
   extern double *intprep_calc;
   FILE *fp;

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
   gets(filename) ; */
   strcpy(dirname+7,filename);
   fp=own_fopen(dirname);
   bldcompos(fp, &compsample);
   fclose(fp);
   convcomp(&compsample, &sample);

   /* now set up and fill block with abscoeff for elemens with labda */

   /* debugging
   printf("Voor call fillmulabda\n");  /* */
   fillmulabda(&sample,  nr_of_meas_lines);
   /* debugging
   printf("Na call fillmulabda\n");  /* */

   /* Now only once for bulk materials */

     calc_alphas(&sample, nr_of_meas_lines);
     /* debug 
     printf("Na calc_alphas\n"); /* */
     conciter(&sample);

     printf("========== after first conciter ======\n");

    /* now simplex to convert to compound composition */

     fitconc(&compsample, &sample);
     convcomp(&compsample, &sample);
     conciter(&sample);

     printf("========== after second conciter ====\n");

     alphrho(&sample, nr_of_meas_lines);
     phase2(&sample);
     fitconc(&compsample, &sample);
     return(0);
}
