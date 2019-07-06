/*
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*	pilabda.c
*
*	calculates primary fluorescence intensity
*	of given element in single  homogeneous bulklayer
*	for monochromatic excitation radiation with 
*       wavelength labda
*
*	formula's from 
*       D.K.G. de Boer et.al. Adv. in X-Ray Anal. 
*	33 (1990) 237
*
*	M.Bos
*	sep 1997
*
*   changed fluorescence condition
*   Jan 2004
*
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
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
extern double edges[TOTAL][10] ;


void show_res(int , double , char *, int ,
              double , double , double ,
              double , double , double ,
              double , double,  double );

double primintfull( double Ci, double Jilambda, double omegai,
                    double gi, double tauilambda, double mu1lambda,
                    double mu1i)
{
	double part1, part2 ;
        part1=Gfactor*Ci*Jilambda*omegai*gi*tauilambda ;
        part2=sin(psi1)*(mu1lambda/sin(psi1)+mu1i/sin(psi2));
        return(part1/part2);
}

double pilabda(int line_nr, int index, layer *player)
{
      double  Ci, Jilambda ;
      double omegai, gi, tauilambda, mu1lambda, mu1i ;
      double  energy, lambdai, lambda, labd0, l_intvl,  Ei;
      double pil, thickness;
      int z, nr_elem_in_layer, hole_nr ;
      layer *ptop;
      char *shellname ;
      char linesymb[6];
      char elsymb[3];


      /* calc lambda belonging to line_nr and index */

     labd0= 12.396/((lines_to_meas[line_nr]).kev);
     l_intvl=(LAMBDA_HIGHEST-labd0)/PNTS_SPECTRUM;
     lambda=labd0+index*l_intvl;
	 energy=12.397/lambda;
     /* printf("In pilbada calctd lambda from index %g\n",
        lambda); /* */

      /* now calculate lamda of the measured line */
      
      strcpy(linesymb,"\0\0\0\0\0\0");
      strcpy(elsymb,"\0\0\0");
      strcpy(linesymb,(lines_to_meas[line_nr]).line);
      strcpy(elsymb, get_asymb((lines_to_meas[line_nr]).z,
                                 &lijst_els));
      lambdai=characwl(linesymb, elsymb,&lijst_els);
/* printf("Na characwl in pilabda for line %s el %s \n",linesymb,elsymb);  /* */
      Ei= 12.396/lambdai;  
 /*     printf(" lambdai %g   \n",lambdai);/* */
 /*     if (lambdai < lambda)  return (0.0); /* removed Jan 2004 */
      z= (lines_to_meas[line_nr]).z ;    /* atoomnummer */
 
      ptop=player;
      thickness= get_mthickness(ptop);
/*      printf("Mass thickn in pilabda found %g\n", thickness); /* */
      nr_elem_in_layer=get_nr_elem(ptop, elsymb);
      Ci=get_weightfract(ptop, nr_elem_in_layer);
      shellname= line_to_shell(linesymb, linenames, shellreeks);
      hole_nr=get_hole_nr(shellname);
	  if (lambda > (12.396/edges[z][hole_nr-1])) return (0.0); /* new cond Jan 2004 */
      Jilambda=getjmpfact(z, hole_nr,energy, jmpdata);
      omegai=get_omega(z, shellname, omegas);
/*      printf("Voor getrelrate z %d line %s \n", z, linesymb);  */
      gi=getrelrate(z, linesymb,relrates);
/*      printf("in pilambda gi van Z %d line %s gi is %e \n",
               z, linesymb, gi); /* */
/*      tauilambda= calc_muilambda(z, energy, &abstabel); /* */
 /*     printf("Line_nr is nu %d  z is %d en index %d\n", line_nr, z, index); 
      /* check_mu_el(z,line_nr); /* */
      tauilambda= muzlambda[line_nr][z][index];
/*       printf("tauilambda is %g \n", tauilambda);
     /* equal to mu no scatt ? */
/*      printf("In pilabda voor mumix_snel\n"); /* */
      mu1lambda=calc_mumix_snel(line_nr, index, ptop, &lijst_els,
                                 &abstabel); 
 /*     printf("In pilabda na mumixsnel \n"); /* */
   
      mu1i=calc_mumix(Ei,ptop, &lijst_els, &abstabel);
/*        show_res(z,Ci,shellname, hole_nr,Jilambda, omegai,gi,tauilambda,
            mu1lambda, lambdai, lambda, Ei, mu1i );/*  */
      pil= primintfull(Ci,Jilambda, omegai, gi, tauilambda,
                       mu1lambda, mu1i);
      pil=pil*selfattn(mu1lambda,mu1i,thickness); 
/*      printf("return reached in pilabda with pil is %g\n", pil); /* */
      return(pil);
}

void show_res(int z, double Ci, char *shellname, int hole_nr,
              double Jilambda, double omegai, double gi,
              double tauilambda, double mu1lambda, double lambdai,
              double lambda, double Ei, double mu1i)
{
  printf("********************************************************\n");
  printf("Z is %d, Ci is %f, shell is %s, hole_nr is %d\n",
          z, Ci,shellname, hole_nr);
  printf("Jilambda is %f, omegai is %f, gi is %f\n",
          Jilambda, omegai, gi);
  printf("tauilambda is %f, mu1lambda is %f, lambdai is %f, lambda is %f\n",
          tauilambda, mu1lambda, lambdai, lambda);
  printf("Ei is %f, mu1i is %f \n",Ei, mu1i);
  printf("********************************************************\n");
}

  

