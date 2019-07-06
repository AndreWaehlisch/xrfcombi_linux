/**********************************************************************
*
*
*    COPYRIGHT (C) M.BOS 1999 , for details see copying
*
*
*                    xrfwang.c
*                    Main program for quantitative 
*                    X-ray fluorescence analysis
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                    
*                    changed to use new lines.dat file
*                    which includes Gi-factor as last entry on line
*
*                    Sept 13, 1999
*
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

#define MAXITERS 1000
int main(int argc,char *argv[])
{
   extern multicomplayer compsample ;   /* sample in terms of compounds */
   extern multilayer sample ;          /* sample in terms of elements */
   extern struct spectrum_pair **pntrspectra ;
   int i, j ;
   extern int ref_line ;
   char el[3];
   char lin[6];
   FILE *fp, *fpconverg, *fpresultaat;
   char filename[80]="sample.dat"; 
   char dirname[80]="./data/" ;
   extern char title[132]; 
   double *free_pars ;  /* pointer to array with simplex params + resp */
   double *pars_min ;
   double *pars_max; 
   int  z_filt;
   double mth_filt, kev;
   int nr_compounds;
   extern double *rik_meas; 
   extern double *rik_calc;
   extern double *intprep_calc;
   extern double *pi_fact;
   extern struct pstruct *p, pcent, **p_p ;
   extern int prt_cycle;
   extern int ndata;
   extern double quad_test, test, ypmin, yzero; 
   extern int quad_cycle, nquad_skip, maxquad_skip;

   /* debugging
   printf("voor init tables \n");  /* */
   init_tables();

   /* debugging 
   printf("na init tables \n");  /* */

   set_spectrom();
   /* debugging 
   printf("na set spectrom \n"); /* */
   fp=own_fopen("./data/lines.dat");
   nr_of_meas_lines=get_nr_lines_meas(fp);
   fseek(fp, 0, SEEK_SET);  /* reset to beginning of lines.dat */

   /* debugging  */
   printf("nr_of_meas_lines %d \n",nr_of_meas_lines); /* */

   rik_meas= malloc(nr_of_meas_lines* sizeof(double));
   rik_calc= malloc(nr_of_meas_lines* sizeof(double));
   intprep_calc= malloc(nr_of_meas_lines* sizeof(double));
   pi_facts= malloc(nr_of_meas_lines* sizeof(double));
   ref_line=get_gifactors(fp); 
   fclose(fp);
   calc_rikmeas(rik_meas, nr_of_meas_lines, ref_line); 
   /* debugging 
   printf("na lezen intensities \n"); /* */
   for (i=0; i < nr_of_meas_lines; i++)
   {
       printf("%d Z: %d line %s Rik  %g \n", i, (lines_to_meas[i]).z,
              (lines_to_meas[i]).line, rik_meas[i]);
   } 


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
   printf("Na call fillmulabda\n"); /* */

   /* Now generate spectra for each line separately */

   for (j=0; j < nr_of_meas_lines; j++)
   {
     strcpy(el, get_asymb((lines_to_meas[j]).z, &lijst_els));
     strcpy(lin, (lines_to_meas[j]).line);
     z_filt=(lines_to_meas[j]).filt_z;
     mth_filt=(lines_to_meas[j]).mth_filt;
     kev=(lines_to_meas[j]).kev; 
     pntrspectra[j]=gen_spec(anode_el, take_off, berryliumd, kev,
                 PNTS_SPECTRUM, z_filt, mth_filt);
     /* debugging 
     printf("Na gen-spec lijn %s el %s\n",lin,el); /* */
   }


/*   printf("Give title for experiment ");
   gets(title); */
   free_pars=malloc( (nr_free_vars+1)*sizeof(double));
   pars_min=malloc( (nr_free_vars+1)*sizeof(double));
   pars_max=malloc( (nr_free_vars+1)*sizeof(double));
   if (pars_max==NULL) 
   { printf("Cannot malloc pars _max \n");
     exit(1);
   }
   /* debugging
   printf("na malloc free_pars \n"); /* */

   /* Now ready for simplex optimisation */

   /* setup first simplex */
   get_params(&compsample,free_pars);
   /* show_params(free_pars);   /*  */

   /* first alloc space needed for the simplex */

   nfree=nr_free_vars ;
   /* debugging
   printf("Nr of free parameters is %d \n",nfree); /* */
   nvert=nfree+1;
   nparm=nfree;
   ndata=nr_of_meas_lines;
   p= (struct pstruct *) malloc( nvert*sizeof(struct pstruct));
   p_p=(struct pstruct ** )malloc( nvert*sizeof (struct pstruct *));
  
   if (p_p==NULL)
   {
	    printf("Cannot malloc p_p \n");
		exit(1);
   }
   /* Now calculate first simplex */

   for (i=0; i < nfree ; i++)
   {
        pars_min[i]=free_pars[i]*0.9;
        pars_max[i]=free_pars[i]*1.1;
   }
 
   for (i=0; i < nvert ; i++)
   {
      p[i].val=0.0;
      for (j=0; j < nfree ; j++)
      {
          p[i].parm[j]=pars_min[j] + ((pars_max[j]-pars_min[j])*i)/nvert;
      }
/*      printf("Voor func call in xrf.c \n"); /* */
      func(&p[i]);
/*      printf("Na func call \n"); /* */
    }

    /* now print start situation with params and ssq */

/*   for (i=0; i < nvert ; i++)
   {
      printf("Vertex %d ",i);
      for (j=0; j < nfree ; j++)
      {
        printf("%g " , p[i].parm[j]);
      }
      printf("ssq %g \n",p[i].val);
    }
/*  */


   exit_test=1.0E-30*ndata ;
   iter=0 ;

   /*  ITERATIONS LOOP ALTERNATING simpfit() and simpdev() */
   
   prt_cycle= 5 ;
   quad_cycle= 5 ;
   nquad_skip = 0;
   quad_test= 5 ;
   maxquad_skip= 25;
   fpconverg=fopen("./data/converge.dat","w");

   while (1)
   {
      maxiter= iter + prt_cycle ;
      simpfit(fpconverg); 
      ffitprint(fpconverg);

      if (iter == maxiter)
      {
         if (quad_test >=1 )
          {
            if (iter < quad_test)
                continue;
          }
         else if ( test > quad_test) 
           continue;
      }
     simpdev(fpconverg);
     fquadprint(fpconverg);

     if (ypmin < yzero)
           nquad_skip=0;
     else if (nquad_skip < maxquad_skip)
          nquad_skip++ ;

     if (quad_test >=1 )
     {
       if (iter >= quad_test)
           quad_test = iter + (nquad_skip + 1)* quad_cycle ;
     }
     else
          if ( test <= quad_test)
              quad_test = quad_test/10 ;
     fflush(fpconverg);
     if ((iter > MAXITERS)|| (test < exit_test)) break ;
    
   }
   fclose(fpconverg);
   if ((fpresultaat=fopen("./data/results.dat","w"))==NULL)
   {
     printf("Cannot open ./data/results\n");
     exit(1);
   }   
     
   if (iter > MAXITERS) printf("Warning maximum # iterations reached \n");

   for (j=0; j < nr_of_meas_lines; j++)
   {
    printf("Line nr %d RXI found %g RXI meas %g\n", j, rik_calc[j],
             rik_meas[j]);
    fprintf(fpresultaat,"Line nr %d RXI found %g RXI meas %g\n", j, rik_calc[j],
             rik_meas[j]);

   }
   nr_compounds=prtmcomp(stdout, &compsample);
   nr_compounds=prtmcomp(fpresultaat, &compsample);
   printf("Ssq is %g\n",pcent.val);
   fprintf(fpresultaat,"Ssq is %g\n",pcent.val);
   fclose(fpresultaat);
   return(0);
}
