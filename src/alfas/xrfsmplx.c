/***********************************************************************

*

*    COPYRIGHT (C) M.BOS 1998 , for details see copying

*

*

*                    xrfsmplx.c

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



#define MAXITERS 400

int main(int argc,char *argv[])

{

   extern multicomplayer compsample ;   /* sample in terms of compounds */

   extern multilayer sample ;          /* sample in terms of elements */

   multilayer bulksample ;

   layer *pbulk;

   double primint_bulk, totint_bulk;

   extern struct spectrum_pair **pntrspectra ;

   double labda_in ;

   int i, j ;

   char el[3];

   char lin[6];

   FILE *fp, *fpconverg;

   char filename[80]="sample.dat"; 

   char dirname[80]="./data/" ;

   extern char title[132]; 

   double *free_pars ;  /* pointer to array with simplex params + resp */

   double *pars_min ;

   double *pars_max; 

   int  z_filt;

   double mth_filt, kev;

   int nr_compounds;

   extern double *rxi_meas; 

   extern double *rxi_calc;

   extern double *intbulk_calc;

   extern double *intprep_calc;

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

   printf("nal lezen rxi's \n"); /* */

/*   for (i=0; i < nr_of_meas_lines; i++)

   {

       printf("%d Z: %d line %s rel int. %g \n", i, (lines_to_meas[i]).z,

              (lines_to_meas[i]).line, rxi_meas[i]);

   }

*/ 





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



   /* Now only once for bulk materials */



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

     pbulk=mkbulklay(el);

     /* now set up complete sample for bulk element */

     bulksample.nr_layers= 1;

     bulksample.lagen=&pbulk; 

     totint_bulk = 0.0;

     for (i=0; i < PNTS_SPECTRUM; i++)

      {

        labda_in=(pntrspectra[j][i]).lambda;

        primint_bulk = pilabda(j, i, pbulk)*(pntrspectra[j][i]).intens;

        totint_bulk += primint_bulk;

      }

    /*  printf("Bulkint %s %s  %e \n",el, lin, totint_bulk); */

     intbulk_calc[j]=totint_bulk;

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

/*   show_params(free_pars);  */



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

*/





   exit_test=2.0E-10*ndata ;

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

 

   if (iter > MAXITERS) printf("Warning maximum # iterations reached \n");

   printf("\n");

   for (j=0; j < nr_of_meas_lines; j++)

   {

    printf("Line  %s %s RXI found %g RXI meas %g\n", 

             get_asymb(lines_to_meas[j].z, &lijst_els),

             lines_to_meas[j].line, rxi_calc[j],

             rxi_meas[j]);

   }

   printf("\n");

   nr_compounds=prtmcomp(stdout, &compsample); 

   printf("Ssq is %g\n",pcent.val);

   return(0);

}

