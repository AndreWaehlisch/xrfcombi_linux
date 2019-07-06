/***********************************************************************
*
*    COPYRIGHT (C) M.BOS 1999 , for details see copying
*
*
*                    fitcomp.c
*                    simplex routine
*                    to fit compound concentrations
*                    on element concentration
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                    September 1999
*
*                    Sep 22 1999 included cndx sum conc =1.0
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

#define MAXITERS 1000

int fitconc(multicomplayer *compsample, multilayer *sample)
{
   int i, j ;
   FILE  *fpconverg;
   extern int nfree, nvert, nparm ;
   extern struct pstruct *p, pcent, **p_p ;
   extern int prt_cycle;
   extern int ndata;
   extern double quad_test, test, ypmin, yzero, exit_test; 
   extern int quad_cycle, nquad_skip, maxquad_skip;
   extern int iter, maxiter, nr_free_vars;
   double *free_pars ;  /* pointer to array with simplex params + resp */
   double *pars_min ;
   double  *pars_max; 

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
   
   get_params(compsample,free_pars);
   /* debugging 
    printf("na get_params \n"); 
    checkcomps(compsample);
    show_params(free_pars);  /* */

   /* first alloc space needed for the simplex */

   nfree=nr_free_vars ;
   /* debugging 
   printf("Nr of free parameters is %d \n",nfree); /* */
   nvert=nfree+1;
   nparm=nfree;
   ndata=nfree+1;
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
      funcnv(&p[i]);
    }

    /* now print start situation with params and ssq 

   for (i=0; i < nvert ; i++)
   {
      printf("Vertex %d ",i);
      for (j=0; j < nfree ; j++)
      {
        printf("%g " , p[i].parm[j]);
      }
      printf("ssq %g \n",p[i].val);
    } /* */



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
      simpfit(fpconverg,funcnv); 
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
     simpdev(fpconverg, funcnv);
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
 
   if (iter > MAXITERS) 
   {
     printf("Warning maximum # iterations reached in fitconc\n"); 
   /*  prtmcomp(stdout, compsample); /*  */ 
     return(1);
   }
   else 
   {
     /*  prtmcomp(stdout, compsample);  /*  */ 
     return(0);
   }
}

