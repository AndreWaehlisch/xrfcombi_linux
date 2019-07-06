/*
*    file raylispfr.c
*
*    routine to calculate 2nd order process
*    for photoelectric interaction followed by 
*    Rayleigh scattering of primary fluorescence radiation
*
*
*    cf G.Tirao and G. Stutz XRS 32 (2003) 13
*
*    M.Bos
*    Nov 29th 2003
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "adapq.h"

extern double Gfactor;
extern double psi1;
extern double psi2;
extern elem_list lijst_els;
extern char *linenames[27];
extern char *shellreeks[27];
extern double jmpdata[ALL_ELMS][2];
extern double omegas[TOTAL][NR_OMEGAS];
extern double relrates[TOTAL][TOT_RELRATES];
extern emudata abstabel;
extern meas_line lines_to_meas[MAX_LINES] ;
extern double muzlambda[MAX_LINES][TOTAL][PNTS_SPECTRUM];
extern double sinpsi1;
extern double sinpsi2;
extern double EEi;
extern layer *ptrlaag;
extern double muci ;
extern double muclabda;
extern double cospsi1 ;
extern double cospsi2;
extern double edges[TOTAL][10];
extern double store_thickness ;

double ispfr(char *elsymb, char *linesymb, double lambda, layer *player, int i, int bulk)
{
  
      double  Ci, Jilambda ;
      double omegai, gi, tauilambda, mu1lambda, mu1i ;
      double  energy, lambdai,  Ei;
      double pil;
      int z, nr_elem_in_layer, hole_nr ;
      layer *ptop;
      char *shellname ;
      double cross_dblint ;
	  static double first_integral_sample ;
      static double first_integral_bulk;
      double part1, part2, part3 ,alfagrens ;
	  double part1a, part2a, part3a ;
      double thickness ;

      energy=12.396/lambda; /* conversion from angstrom to eV */

      /* now calculate lamda of the measured line */
      lambdai=characwl(linesymb, elsymb,&lijst_els);
/* printf("Na characwl in ispfr for line %s el %s \n",linesymb,elsymb);  /* */
      Ei= 12.396/lambdai;  
/*      printf(" lambdai %g  lambda %g \n",lambdai, lambda); /* */
     /* if (lambdai < lambda)  return (0.0); /* line is not fluoresced  not correct ?*/
      z= get_z(elsymb, &lijst_els) ;    /* atoomnummer */
      ptop=player;
      thickness= get_mthickness(ptop);
      nr_elem_in_layer=get_nr_elem(ptop, elsymb);
      Ci=get_weightfract(ptop, nr_elem_in_layer);
      shellname= line_to_shell(linesymb, linenames, shellreeks);
      hole_nr=get_hole_nr(shellname);
   	 if (lambda > (12.396/edges[z][hole_nr-1])) return (0.0); /* not fluoresced */
      Jilambda=getjmpfact(z, hole_nr,energy, jmpdata);
      omegai=get_omega(z, shellname, omegas);
/*      printf("Voor getrelrate z %d line %s \n", z, linesymb);  */
      gi=getrelrate(z, linesymb,relrates);
/*      printf("in ispfr gi van Z %d line %s gi is %e \n",
               z, linesymb, gi);  */
      tauilambda= calc_muilambda(z, energy, &abstabel);

     /* equal to mu no scatt ? */
       mu1lambda=calc_mumix(energy, ptop, &lijst_els, &abstabel); 

   
      mu1i=calc_mumix(Ei,ptop, &lijst_els, &abstabel);
/*       show_res(z,Ci,shellname, hole_nr,Jilambda, omegai,gi,tauilambda,
            mu1lambda, lambdai, lambda, Ei, mu1i );  /* */
      pil= primintfull(Ci,Jilambda, omegai, gi, tauilambda,
                       mu1lambda, mu1i);
      /* now set globals for parameters for easy access
         in double integral routine  */

      EEi=Ei;
	  sinpsi1=sin(psi1);
      sinpsi2=sin(psi2);
      ptrlaag=player;
      muci=mu1i;
      muclabda=mu1lambda;
/*	  printf("muclabda %g\n",muclabda); /*  */
      cospsi1=cos(psi1);
      cospsi2=cos(psi2);
      store_thickness = thickness;

	  alfagrens = muci*sinpsi1/muclabda;
/*	  printf("muclabda is %g\n",muclabda);  /* */

	  if (i==0)
	  {
        if (alfagrens >= 1 )
		{
		    alfagrens = PI/2 ;	      
			part1=aquad2d(0.0, alfagrens-1.0e-10, 0.0, PI*2, 1e-6, (ext2fun) fotrayl_func1,
		         (void*) NULL) ;  /*  methode Jan */

		}
		else
		{
		   alfagrens = acos(alfagrens);
           part1=aquad2d(0.0, alfagrens-1.0e-10, 0.0, PI*2, 1e-6,
				(ext2fun) fotrayl_func1,(void*) NULL) ;  /*  methode Jan */
           part1a=aquad2d(alfagrens+1.0e-10, PI/2, 0.0, PI*2, 1e-6, (ext2fun) fotrayl_func1,
		         (void*) NULL) ;  /*  methode Jan */
/*		   printf("part1 is %g part 1a is %g \n", part1, part1a);  /*  */
		   part1 = part1 + part1a ;
		}

		  if (bulk==1)
		  {
			  first_integral_bulk=part1;
		  }
		  else
		  {
			  first_integral_sample=part1;
		  }
	  }
	  else
	  {
		  
		  if (bulk==1)
		  {
			  part1=first_integral_bulk;
		  }
		  else
		  {
			  part1=first_integral_sample;
		  }
	  }
	   alfagrens= acos(-sinpsi2);
/*	   printf("alfagrens is %g \n",alfagrens) ; /* */
       part2 =aquad2d(PI, alfagrens+1.0e-10, 0.0, 2*PI, 1e-6, (ext2fun) fotrayl_func2,
					 (void*) NULL);  /* van Jan !! */
/*	   printf("with pi to grens part2 is %g \n", part2 ); /* */
	   part3 = aquad2d(alfagrens-1.0e-10, PI/2+1.0e-10,  0.0, 2*PI, 1e-6, (ext2fun) fotrayl_func2,
		              (void*) NULL) ; 
/*	   printf("with grens to 0 %g\n", part3); /* */
       cross_dblint= part1+part2+part3 ;
/*	   printf("crossdbl %g \n", cross_dblint);  /* */
     /*  cross_dblint=part1+simpsdbl(integr_func2, PI, PI/2, 0.0, PI*2,NMAX1,NMAX2); /* */

    
     cross_dblint *=pil;

/*   printf("ispfr lambda %g with part1 %g cros_dblint is %g\n", lambda, part1, cross_dblint); /* */
      return(cross_dblint);
}
