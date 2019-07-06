/*
*    file raylisesr.c
*
*    routine to calculate 2nd order process
*    for rayleigh interaction followed by 
*    fotoelectric interaction of scattered radiation
*
*    Now cf S.  Mori and M.Mantler Adv. Xray Anal. 36 (1993) 47
*
*    M.Bos
*    Jan 12th 2004
*
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "adapq.h"
#include "xrfluor.h"

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
extern  double store_thickness ;

double raylisesr(char *elsymb, char *linesymb, double lambda, layer *player, int i, int bulk)
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
      double part1, part2, alfagrens, alfagrens1;

      energy=12.396/lambda; /* conversion from angstrom to eV */
	  
      /* now calculate lamda of the measured line */
      lambdai=characwl(linesymb, elsymb,&lijst_els);
     /* printf("Na characwl in ispfr for line %s el %s \n",linesymb,elsymb);  /* */
      Ei= 12.396/lambdai;  
     /* printf(" lambdai %g  lambda %g \n",lambdai, lambda); /* */
      z= get_z(elsymb, &lijst_els) ;    /* atoomnummer */
      ptop=player;
      thickness= get_mthickness(ptop);
	  store_thickness = thickness;
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
      pil= Jilambda*omegai*gi*Ci*tauilambda ;
	  pil /= sin(psi1)*(mu1lambda/sin(psi1)+mu1i/sin(psi2));

      /* now set globals for parameters for easy access
         in double integral routine  */
      EEi=Ei;
    /*  printf("EEi is %g \n", EEi);  /* */
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
	/*  printf("store_energy is %g muclabda %g muci %g \n",store_energy, muclabda, muci); /* */

      alfagrens = sinpsi1 ;

	  part1 = aquad2d(0.0, alfagrens-1.0e-10, 0.0 , PI*2, 1e-6, (ext2fun) raylfot_func1,
		         (void*) NULL) ;  /*  methode Jan */
      part1 += aquad2d(alfagrens+1.0e-10, 1.0, 0.0, PI*2, 1e-6, (ext2fun) raylfot_func1,
		         (void*) NULL) ;  /*  methode Jan */

	  alfagrens1 = -(muclabda/muci)*sinpsi2 ;
	  if (alfagrens1 < -1)
	  {
         part2 =aquad2d(-1.0, -1.0e-10, 0.0, 2*PI, 1e-6, (ext2fun) raylfot_func2,
					 (void*) NULL);  /* van Jan !! */
	  }
	  else
	  {
		  alfagrens = alfagrens1;
	/*	  printf("alfagrens %g\n", alfagrens);   /*  */
		  part2 = aquad2d(-1.0 , alfagrens-1.0e-10 , 0.0, 2*PI, 1e-6, (ext2fun) raylfot_func2,
					 (void*) NULL);  /* van Jan !! */
          part2 += aquad2d(alfagrens+1.0e-10, -1.0e-10, 0.0, 2*PI, 1e-6, (ext2fun) raylfot_func2,
					 (void*) NULL);  /* van Jan !! */
	  }

       cross_dblint = part1 + part2 ;
	/*   printf("part1 %g part2 %g crossdbl %g \n", part1,part2,  cross_dblint);   /*  */
    
     cross_dblint *= pil;

 /*  printf("raylisesr lambda %g with part1 %g cros_dblint is %g\n", lambda, part1, cross_dblint); /* */
      return(cross_dblint);
}
