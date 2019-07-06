/*
*       COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*	secfluijscat.c
*
*	calculates secundary fluorescence intensity
*	of given element in single  homogeneous bulklayer
*	layer is radiated with monochromatic  radiation with 
*       wavelength labda and the secundary fluorescence
*       caused by the radiation of line lambda_j 
*       of element j is calculated
*
*
*	formula's from 
*       D.K.G. de Boer et.al. Adv. in X-Ray Anal. 
*	33 (1990) 237
*
*	M.Bos
*	sep 1997
*
*/
#include <stdio.h>
#include <math.h>
#include "xrfluor.h"

extern double psi1;
extern double psi2;
extern elem_list lijst_els;
extern char *linenames[27];
extern char *shellreeks[27];
extern double jmpdata[ALL_ELMS][2];
extern double omegas[TOTAL][NR_OMEGAS];
extern double relrates[TOTAL][TOT_RELRATES];
extern double kev;
extern emudata abstabel;



void show_res_scat(int , double , char *, int ,
              double , double , double ,
              double , double , double ,
              double , double, double );

double secijlambda_scat( double Cj, double Jjlambda, double omegaj,
                    double gj, double Jilambdaj, double Jilambda,
                    double taujlambda, double tauilambda,
                    double tauij, double mu1i, double mu1lambda,
                    double mu1j, double thickness, double Pilambda)
{
	double part1, part2 , part3, mu1iprime, mu1lambdaprime ;
        part1= 0.5*Cj*Jjlambda*omegaj*gj ;
        part2= (Jilambdaj/Jilambda)*((taujlambda*tauij)/(tauilambda*mu1j));
        mu1iprime= mu1i/sin(psi2);
        mu1lambdaprime= mu1lambda/sin(psi1);
        part3= ltot(mu1iprime,mu1lambdaprime,mu1j,thickness)*Pilambda; 
/*        part3= linfty(mu1iprime,mu1lambdaprime,mu1j)*Pilambda; */
        return(part1*part2*part3);
}



double sijlambda_scat(char *el_i, char *el_j, char *line_i, char *line_j,
                 layer *ptop, double lambda, double Pilambda)
{
   double Cj, Jjlambda, omegaj, gj, Jilambdaj, Jilambda;
   double taujlambda, tauilambda, tauij, mu1i, mu1lambda ;
   double mu1j, thickness,  energy_in, lambdaj ;
   double Ej, lambdai,  Ei,    result;    
   int zi , zj, nr_elem_enhancing, hole_nr ;
   char *shellname_i, *shellname_j;

   energy_in=12.396/lambda;
   zi= get_z(el_i, &lijst_els);   /* atom number measured element */
   zj= get_z(el_j, &lijst_els);   /* enhancing element */
   thickness=get_mthickness(ptop);
   nr_elem_enhancing= get_nr_elem(ptop, el_j);
   Cj=get_weightfract(ptop, nr_elem_enhancing);

   shellname_i= line_to_shell(line_i, linenames, shellreeks);
   hole_nr=get_hole_nr(shellname_i);
   Jilambda=getjmpfact(zi, hole_nr, energy_in, jmpdata);
   /* printf(" Jilambda meas = %g for Z = %d, hole %d energy_in %g\n",
                 Jilambda, zi, hole_nr, energy_in); */

   shellname_j= line_to_shell(line_j, linenames, shellreeks);
   hole_nr=get_hole_nr(shellname_j);
   Jjlambda=getjmpfact(zj, hole_nr, energy_in, jmpdata);
   /* printf("Jjlambda enhanc = %g for Z = %d, hole %d energy_in %g\n",
             Jjlambda, zj, hole_nr, energy_in); */

   hole_nr=get_hole_nr(shellname_i);
   lambdaj=characwl(line_j, el_j, &lijst_els); 
   Ej= 12.396/lambdaj;
   Jilambdaj=getjmpfact(zi, hole_nr, Ej, jmpdata);  
   /* printf("Jilambdaj = %g for Z %d hole_nr %d energy_in %g \n",
            Jilambdaj, zi, hole_nr, Ej);  */
   taujlambda= calc_muilambda(zj, energy_in, &abstabel);
   tauilambda= calc_muilambda(zi, energy_in, &abstabel);

   /* now calc energ of fluorescent radiation of element j */
   
   lambdaj=characwl(line_j, el_j, &lijst_els);
   Ej=12.396/lambdaj ;
   tauij= calc_muilambda(zi, Ej, &abstabel);
 
   lambdai=characwl(line_i, el_i, &lijst_els);
   Ei=12.396/lambdai ;

   mu1i= calc_mumix(Ei, ptop, &lijst_els, &abstabel);
   mu1j= calc_mumix(Ej, ptop, &lijst_els, &abstabel);
   mu1lambda = calc_mumix(energy_in, ptop, &lijst_els, &abstabel);
   omegaj=get_omega(zj, shellname_j, omegas);
   gj=getrelrate(zj, line_j, relrates);
   result=secijlambda_scat(Cj, Jjlambda, omegaj, gj, Jilambdaj, Jilambda,
                      taujlambda, tauilambda, tauij, mu1i, mu1lambda,
                      mu1j, thickness, Pilambda);

/************************************
*   debugging 
   printf("Cj is %f \n", Cj);
   printf("Jjlambda is %e \n", Jjlambda);
   printf("omegaj is %e \n",omegaj);
   printf("gj is %e \n",gj);
   printf("Jilambda is %e \n", Jilambda);
   printf("Jilambdaj is %e \n",Jilambdaj);
   printf("taujlambda is %e \n", taujlambda);
   printf("tauilambda is %e \n", tauilambda);
   printf("tauij os %e \n",tauij);
   printf("mu1i is %e \n", mu1i);
   printf("mu1lambda is %e \n", mu1lambda);
   printf("mu1j is %e \n", mu1j);
   printf("thickness is %e \n", thickness);
   printf("Pilambda is %e \n", Pilambda);
   */ 
return(result);
}

