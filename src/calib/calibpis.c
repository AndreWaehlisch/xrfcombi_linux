/***********************************************************************
*
*   copyright (c) M.Bos, for details see copying
*
*
*                    calibpis.c
*                    Program for quantitative 
*                    X-ray fluorescence analysis
*                    calibration of measured intensities
*                    versus total intensities (for method Wang)
*                    based on fundamental parameter approach
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                    Januray 1999
*
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

extern FILE *own_fopen(char *);

int  main(int argc, char *argv[])
{
   multilayer bulksample ;
   layer *ptop, *pbulk;
   double primint, secint, totint, primint_bulk, totint_bulk;
   struct spectrum_pair *pntrspectrum ;
   double labda_in ;
   int i, j, k ;
   char el[3];
   char lin[6];
   int z_filt;
   double mth_filt, kev;
   FILE *fpsmpl, *fpcalres ;
   int nr_smpls;
   multilayer samples[MAX_SMPLS];
   char filenaam[100]="./data/concset.dat";
   multicomplayer monsters[MAX_SMPLS];
   double **rxi ;
   double **meas_ints;
   double *sigmxy;
   double *sigmxx ;
   double *calib_slope, *sigma2_slope;
   
   double kmin; 
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
   }  /* */

   /* printf("Give name of file with sample data ");  
   gets(filenaam); */

   fpsmpl=own_fopen(filenaam);


   fscanf(fpsmpl, "%d", &nr_smpls);
   rxi= (double **) malloc(nr_smpls*(sizeof( double *)));
   meas_ints= (double **)malloc(nr_smpls*(sizeof(double *)));
   /* debugging
   printf("nr_smpls is %d \n", nr_smpls); */
   for (k=0; k < nr_smpls; k++)
   {
     /* debugging
     printf("K-loop k is % d\n",k); /* */
     rxi[k]=(double *)malloc(nr_of_meas_lines*(sizeof(double)));
     meas_ints[k]=(double *)malloc(nr_of_meas_lines*(sizeof(double)));
     if (meas_ints[k]==NULL)
     {
       printf("error in malloc \n");
       exit(1);
     }
     def_sample(fpsmpl, &samples[k], &monsters[k]);
     ptop=get_layer(&samples[k], 0);
     for (j=0; j < nr_of_meas_lines; j++)
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
       pbulk=mkbulklay(el);
       /* now set up complete sample for bulk element */
       bulksample.nr_layers= 1;
       bulksample.lagen=&pbulk; 
       totint=0.0;
       totint_bulk = 0.0;
       for (i=0; i < PNTS_SPECTRUM; i++)
       {
        labda_in=(pntrspectrum[i]).lambda;
/* printf("Voor pilamba i is %d en intens %g \n",i, (pntrspectrum[i]).intens);  /* debugging */
        primint = pilabda(el,lin,labda_in, ptop)*(pntrspectrum[i]).intens; 
        primint_bulk = pilabda(el,lin,labda_in, pbulk)*(pntrspectrum[i]).intens;
        secint = intrasecflu(el,lin, labda_in , &samples[k], 0, primint);
    /*    printf("primint %g primbulk %g secint %g\n",
                primint, primint_bulk, secint);  /* */
      /*  printf("secint %g \n",secint); /* */
        totint += primint+secint;
        totint_bulk += primint_bulk;
     /*   printf(" totint %g totint_bulk %g \n", totint, totint_bulk); /* */
       }
       /*
       printf("totint %s %s is %e bulkint %e\n",el, lin, totint,totint_bulk);
       printf("Rel. int %s %s is %e \n",el,lin, totint/totint_bulk); /* */
       rxi[k][j]=totint; /* here for rxi divide by totint_bulk */
     }
    }
    if (fclose(fpsmpl)== EOF) 
    {
      printf("Problems with closing file \n");
      exit(1);
    } ;

  /* get measured intensities */
  /*   printf("Give name of file with measured intensities ");
   gets(filenaam); */
   strcpy(filenaam, "./data/intset.dat");
   fpsmpl=own_fopen(filenaam);
   fpcalres=fopen("./data/calibwang.dat","w");
   if (fpcalres==NULL)
   {
     printf("Cannot open file ./data/calibres.dat\n");
     exit(1);
   }

   
   /* debugging
   printf("voor lezen meas ints nr_smpls is %d \n", nr_smpls);  /* */
   for (i=0; i < nr_smpls; i++)
   {
     for (j=0; j < nr_of_meas_lines; j++)
     fscanf(fpsmpl,"%lf",&meas_ints[i][j]);
   }
   /* debugging  
  printf("Na lezen intensiteiten nr_smpls is %d\n", nr_smpls);  /* */

   sigmxy= (double *)malloc(nr_of_meas_lines*(sizeof(double)));
   sigmxx= (double *)malloc(nr_of_meas_lines*(sizeof(double)));
   calib_slope = (double *)malloc(nr_of_meas_lines*(sizeof(double)));
   sigma2_slope = (double *)malloc(nr_of_meas_lines*(sizeof(double)));

   /* calc calibration line  */
   for (j=0; j < nr_of_meas_lines; j++)
   {
    sigmxy[j]=0.0;
    sigmxx[j]=0.0;
    for (i=0; i < nr_smpls ; i++)
    {
     /* debugging */
     printf("sample %d rxi %g meas %g\n", i, rxi[i][j],meas_ints[i][j]); /* */
      sigmxy[j] += rxi[i][j]*meas_ints[i][j];
      sigmxx[j] += rxi[i][j]*rxi[i][j];
    }
    calib_slope[j]=sigmxy[j]/sigmxx[j] ;
    fprintf(fpcalres,"%e\n",calib_slope[j]);
   }
   fflush(fpcalres);
   printf("Na flush\n");
/*   if (fclose(fpcalres)== EOF)
   {
      printf("Problems with closing file \n");
      exit(1);
   } */
   printf("na close fpcalres \n");
   fclose(fpsmpl);
   /* debugging */
   printf("Na printen calibres.dat \n");


   fpcalres=fopen("./data/stdevcal.dat","w");
   if (fpcalres==NULL)
   {
     printf("Cannot open file ./data/stdevcal.dat\n");
     exit(1);
   }
   for (j=0; j < nr_of_meas_lines; j++)
   {
     kmin=0.0;
     for (i=0; i < nr_smpls; i++)
     {
        kmin+= (calib_slope[j]*rxi[i][j]-meas_ints[i][j])*
                (calib_slope[j]*rxi[i][j]-meas_ints[i][j]);
     }
     sigma2_slope[j]= (kmin/(nr_smpls-1))/sigmxx[j];
     fprintf(fpcalres,"Lijn nr %d stdev slope %e procent\n", j, 
              100*sqrt(sigma2_slope[j])/calib_slope[j]);
   }  
   fflush(fpcalres);
   if (fclose(fpcalres)== EOF)
   {
      printf("Problems with closing file \n");
      exit(1);
   } 

  /* print results for gnuplot */


   for (j=0; j < nr_of_meas_lines ; j++)
   {
     printf("%g \n",calib_slope[j]);
     sprintf(filenaam,"./data/calib%1d.dat",j);
     fpsmpl=fopen(filenaam,"w");
     for (i=0; i < nr_smpls ; i++)
     {
       fprintf(fpsmpl,"%g %g %g\n",rxi[i][j], meas_ints[i][j], rxi[i][j]*calib_slope[j]);
     }
     fclose(fpsmpl);
   }     
  return 0;
}
