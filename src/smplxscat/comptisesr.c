/*
*    file comptisesr.c
*
*    routine to calculate 2nd order process
*    for compton interaction followed by 
*    fotoelectric interaction of scattered radiation
*
*    cf G. Tirao and G. Stutz
*    XRS 32(2003)13 
*
*    M.Bos
*    Jan 17th 2004
*
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
extern double store_labda ;
extern double store_labdai ;
extern int store_z ;
extern double store_tauilambda ;
extern double store_energy ;
extern double store_thickness ;
extern double fluor_grens;

double comptisesr(char *elsymb, char *linesymb, double lambda, layer *player, int i, int bulk)
{
  
      double  Ci, Jilambda ;
      double omegai, gi, tauilambda, mu1lambda, mu1i ;
      double  energy, lambdai,  Ei;
      double pil, thickness;
      int z, nr_elem_in_layer, hole_nr ;
      layer *ptop;
      char *shellname ;
      double cross_dblint ;
	  static double first_integral_sample ;
      static double first_integral_bulk;
      double part1, alfagrens ;

      energy=12.396/lambda; /* conversion from angstrom to eV */
	  
      /* now calculate lamda of the measured line */
      lambdai=characwl(linesymb, elsymb,&lijst_els);
/*      printf("Na characwl in ispfr for line %s el %s \n",linesymb,elsymb);  /* */
      Ei= 12.396/lambdai;  
  /*    printf(" lambdai %g  lambda %g \n",lambdai, lambda); /* */
      z= get_z(elsymb, &lijst_els) ;    /* atoomnummer */
      ptop=player;
      thickness= get_mthickness(ptop);
      nr_elem_in_layer=get_nr_elem(ptop, elsymb);
      Ci=get_weightfract(ptop, nr_elem_in_layer);
      shellname= line_to_shell(linesymb, linenames, shellreeks);
      hole_nr=get_hole_nr(shellname);
	  fluor_grens=12.396/edges[z][hole_nr-1];
	  if (lambda > fluor_grens ) return (0.0); /* not fluoresced */
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
/*	  show_res_scat(z,Ci,shellname, hole_nr,Jilambda, omegai,gi,tauilambda,
            mu1lambda, lambdai, lambda, Ei, mu1i );  /* */
      pil= Jilambda*omegai*gi*Ci ;
	  pil /= sin(psi1)*(mu1lambda/sin(psi1)+mu1i/sin(psi2));

      /* now set globals for parameters for easy access
         in double integral routine  */
      EEi=Ei;
/*      printf("EEi is %g \n", EEi);  /* */
	  sinpsi1=sin(psi1);
      sinpsi2=sin(psi2);
      ptrlaag=player;
      muci=mu1i;
      muclabda=mu1lambda;
      cospsi1=cos(psi1);
      cospsi2=cos(psi2);
	  store_labda = lambda ;
	  store_labdai = lambdai ;
	  store_z = z;
      store_energy = energy;
	  store_thickness = thickness;
/*	  printf("store_energy is %g\n",store_energy); /* */

      
	  part1=aquad2d(0.0+1.0e-10, 1.0 , 0.0, PI*2, 1e-6, (ext2fun) comptfot_func1,
		         (void*) NULL) ;  /*  methode Jan */
       cross_dblint=part1+aquad2d(-1.0 , 0.0, 0.0, 2*PI, 1e-6, (ext2fun) comptfot_func2,
					 (void*) NULL);  /* van Jan !! */

     cross_dblint *= pil;

 /* printf("comptisesr lambda %g with part1 %g cros_dblint is %g\n", lambda, part1, cross_dblint); /* */
      return(cross_dblint);
}

