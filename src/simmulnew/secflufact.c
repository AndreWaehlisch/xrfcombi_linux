/*
*       COPYRIGHT M.BOS (c) 1998, for details see file copying
*
*	secflufact.c
*
*	calculates factor in secundary fluorescence intensity
*       equation
*
*	formula's from 
*       D.K.G. de Boer et.al. Adv. in X-Ray Anal. 
*	33 (1990) 237
*
*	M.Bos
*	sep 1998
*
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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




double secflufact( double Cj, double Jjlambda, double omegaj,
                    double gj, double Jilambdaj, double Jilambda,
                    double taujlambda, double tauilambda,
                    double tauij)
{
	double part1, part2  ;
        part1= 0.5*Cj*Jjlambda*omegaj*gj ;
        part2= (Jilambdaj/Jilambda)*((taujlambda*tauij)/(tauilambda));
        /* debugging 
        printf("In secflufact part1 %g part2 %g \n", part1, part2); /* */
        return(part1*part2);
}


double sij_lambda_up(char *el_i, char *el_j, char *line_i, char *line_j,
                 multilayer *sample,
                 int nr_curr_layer, int nr_other_layer, double lambda,
                  double primint)
{
   double Cj, Jjlambda, omegaj, gj, Jilambdaj, Jilambda;
   double taujlambda, tauilambda, tauij, mu1i, mu1lambda ;
   double mu1j, energy_in, lambdaj ;
   double thickn_b, mubj, mubl ;
   double Ej, lambdai,  Ei,    result;    
   int zi , zj, nr_elem_enhancing, hole_nr ;
   char *shellname_i, *shellname_j;
   layer *pcurrlayer, *potherlayer ;
   double d1, d2 ;
   double mu2l,  mu2j ;
   double tmpnep ;


   pcurrlayer = sample->lagen[nr_curr_layer];
   potherlayer= sample->lagen[nr_other_layer];
   nr_elem_enhancing = get_nr_elem(potherlayer, el_j);
   
   /* debugging 
   printf("enhancing element is nr %d \n", nr_elem_enhancing); /* */

   /* debugging 
   checkdata(sample);  /* */

   Cj= get_weightfract(potherlayer, nr_elem_enhancing);

   d1= (*pcurrlayer).massthickness ;
   d2= (*potherlayer).massthickness ;

   energy_in=12.396/lambda;
   zi= get_z(el_i, &lijst_els);   /* atom number measured element */
   zj= get_z(el_j, &lijst_els);   /* enhancing element */

   shellname_i= line_to_shell(line_i, linenames, shellreeks);
   hole_nr=get_hole_nr(shellname_i);
   Jilambda=getjmpfact(zi, hole_nr, energy_in, jmpdata);
   /* debugging 
   printf(" Jilambda meas = %g for Z = %d, hole %d energy_in %g\n",
                 Jilambda, zi, hole_nr, energy_in);  /* */

   shellname_j= line_to_shell(line_j, linenames, shellreeks);
   hole_nr=get_hole_nr(shellname_j);
   Jjlambda=getjmpfact(zj, hole_nr, energy_in, jmpdata);
   /* debugging 
   printf("Jjlambda enhanc = %g for Z = %d, hole %d energy_in %g\n",
             Jjlambda, zj, hole_nr, energy_in); /*  */

   hole_nr=get_hole_nr(shellname_i);
   lambdaj=characwl(line_j, el_j, &lijst_els); 
   Ej= 12.396/lambdaj;
   Jilambdaj=getjmpfact(zi, hole_nr, Ej, jmpdata);  
   /* debugging 
   printf("Jilambdaj = %g for Z %d hole_nr %d energy_in %g \n",
            Jilambdaj, zi, hole_nr, Ej);  /* */
   taujlambda= calc_muilambda(zj, energy_in, &abstabel);
   tauilambda= calc_muilambda(zi, energy_in, &abstabel);

   /* now calc energ of fluorescent radiation of element j */
   
   lambdaj=characwl(line_j, el_j, &lijst_els);
   Ej=12.396/lambdaj ;
   tauij= calc_muilambda(zi, Ej, &abstabel);
 
   lambdai=characwl(line_i, el_i, &lijst_els);
   Ei=12.396/lambdai ;
   if (abs(nr_curr_layer - nr_other_layer) > 1)
   {
     thickn_b = calc_th_b(sample, nr_curr_layer, nr_other_layer);
     mubl= calc_mean_mu(sample, nr_curr_layer, nr_other_layer, lambda);
     mubj= calc_mean_mu(sample, nr_curr_layer, nr_other_layer,lambdaj);
     /* debugging 
     printf("curr %d other %d lambdaj %g mubj %g\n", nr_curr_layer,
             nr_other_layer, lambdaj, mubj); /* */
   }
   else
   {
     thickn_b = 0.0;
     mubl=0.0 ;
     mubj =0.0;
   }
   mu2l= calc_mumix(energy_in, potherlayer, &lijst_els, &abstabel);
   mu2j= calc_mumix(Ej, potherlayer, &lijst_els, &abstabel);
   mu1j= calc_mumix(Ej, pcurrlayer, &lijst_els, &abstabel);


   mu1i= calc_mumix(Ei, pcurrlayer, &lijst_els, &abstabel);
   mu1lambda = calc_mumix(energy_in, pcurrlayer, &lijst_els, &abstabel);
   omegaj=get_omega(zj, shellname_j, omegas);
   gj=getrelrate(zj, line_j, relrates);


   result=secflufact(Cj, Jjlambda, omegaj, gj, Jilambdaj, Jilambda,
                      taujlambda, tauilambda, tauij);
   /* debugging 
   printf(" secflufact is %g \n", result);  /* */
   tmpnep = mup(mu1lambda/sin(psi1), mu1i/sin(psi2), mubl/sin(psi1),
                mu2l/sin(psi1),  mu1j, mu2j, mubj, d1, d2, thickn_b) ;

   /* debugging 
   printf("returnwaarde mup is %g \n", tmpnep);  /* */

         result *= tmpnep ;
 
   /* debugging 
   printf("result van sij_lambda_up is %g \n", result);    /* */
    

/************************************
*   debugging  
   printf("========================\n");
   printf("Cj is %f \n", Cj);
   printf("Jjlambda is %e \n", Jjlambda);
   printf("omegaj is %e \n",omegaj);
   printf("gj is %e \n",gj);
   printf("Jilambda is %e \n", Jilambda);
   printf("Jilambdaj is %e \n",Jilambdaj);
   printf("taujlambda is %e \n", taujlambda);
   printf("tauilambda is %e \n", tauilambda);
   printf("tauij os %e \n",tauij);
   printf("d1 is %g \n", d1);
   printf("d2 is %g \n", d2);
   printf("thickn_b is %g \n", thickn_b);
   printf("zi is %d \n", zi);
   printf("zj is %d \n", zj);
   printf("lambdaj is %g \n", lambdaj);
   printf("lambdai is %g \n", lambdai);
   printf("mubl %g mubj %g \n", mubl, mubj);
   printf("mu2l %g mu2j %g mu1j %g mu1i %g mu1lambda %g",
           mu2l, mu2j, mu1j, mu1i, mu1lambda);
   printf("===================\n");
    /*  */ 
return(result);
}


double sij_lambda_down(char *el_i, char *el_j, char *line_i, char *line_j,
                 multilayer *sample,
                 int nr_curr_layer, int nr_other_layer, double lambda,
                  double primint)
{
   double Cj, Jjlambda, omegaj, gj, Jilambdaj, Jilambda;
   double taujlambda, tauilambda, tauij, mu1i, mu1lambda ;
   double mu1j,  energy_in, lambdaj ;
   double thickn_b, mubj, mubl ;
   double Ej, lambdai,  Ei,    result;    
   int zi , zj, nr_elem_enhancing, hole_nr ;
   char *shellname_i, *shellname_j;
   layer *pcurrlayer, *potherlayer ;
   double d1, d2 ;
   double mu2l,  mu2j ;
   double tmpnep ;


   pcurrlayer = sample->lagen[nr_curr_layer];
   potherlayer= sample->lagen[nr_other_layer];
   nr_elem_enhancing = get_nr_elem(potherlayer, el_j);

   /* debugging 
   printf("enhancing element is nr %d \n", nr_elem_enhancing);  /* */

   /* debugging 
   checkdata(sample);  /* */

   Cj= get_weightfract(potherlayer, nr_elem_enhancing);

   d1= (*pcurrlayer).massthickness ;
   d2= (*potherlayer).massthickness ;

   energy_in=12.396/lambda;
   zi= get_z(el_i, &lijst_els);   /* atom number measured element */
   zj= get_z(el_j, &lijst_els);   /* enhancing element */

   shellname_i= line_to_shell(line_i, linenames, shellreeks);
   hole_nr=get_hole_nr(shellname_i);
   Jilambda=getjmpfact(zi, hole_nr, energy_in, jmpdata);
   /* debugging 
   printf(" Jilambda meas = %g for Z = %d, hole %d energy_in %g\n",
                 Jilambda, zi, hole_nr, energy_in);  /* */

   shellname_j= line_to_shell(line_j, linenames, shellreeks);
   hole_nr=get_hole_nr(shellname_j);
   Jjlambda=getjmpfact(zj, hole_nr, energy_in, jmpdata);
   /* debugging 
   printf("Jjlambda enhanc = %g for Z = %d, hole %d energy_in %g\n",
             Jjlambda, zj, hole_nr, energy_in); /*  */

   hole_nr=get_hole_nr(shellname_i);
   lambdaj=characwl(line_j, el_j, &lijst_els); 
   Ej= 12.396/lambdaj;
   Jilambdaj=getjmpfact(zi, hole_nr, Ej, jmpdata);  
   /* debugging 
   printf("Jilambdaj = %g for Z %d hole_nr %d energy_in %g \n",
            Jilambdaj, zi, hole_nr, Ej);  /* */
   taujlambda= calc_muilambda(zj, energy_in, &abstabel);
   tauilambda= calc_muilambda(zi, energy_in, &abstabel);

   /* now calc energ of fluorescent radiation of element j */
   
   lambdaj=characwl(line_j, el_j, &lijst_els);
   Ej=12.396/lambdaj ;
   tauij= calc_muilambda(zi, Ej, &abstabel);
 
   lambdai=characwl(line_i, el_i, &lijst_els);
   Ei=12.396/lambdai ;
   if (abs(nr_curr_layer - nr_other_layer) > 1)
   {
     thickn_b = calc_th_b(sample, nr_curr_layer, nr_other_layer);
     mubl= calc_mean_mu(sample, nr_curr_layer, nr_other_layer, lambda);
     mubj= calc_mean_mu(sample, nr_curr_layer, nr_other_layer,lambdaj);
   }
   else
   {
     thickn_b = 0.0;
     mubl=0.0 ;
     mubj =0.0;
   }
   mu2l= calc_mumix(energy_in, potherlayer, &lijst_els, &abstabel);
   mu2j= calc_mumix(Ej, potherlayer, &lijst_els, &abstabel);
   mu1j= calc_mumix(Ej, pcurrlayer, &lijst_els, &abstabel);


   mu1i= calc_mumix(Ei, pcurrlayer, &lijst_els, &abstabel);
   mu1lambda = calc_mumix(energy_in, pcurrlayer, &lijst_els, &abstabel);
   omegaj=get_omega(zj, shellname_j, omegas);
   gj=getrelrate(zj, line_j, relrates);




   result=secflufact(Cj, Jjlambda, omegaj, gj, Jilambdaj, Jilambda,
                      taujlambda, tauilambda, tauij);
   /* debugging 
   printf(" secflufact is %g \n", result);  /* */
   tmpnep = mdown(mu1lambda/sin(psi1), mu1i/sin(psi2), mubl/sin(psi1),
                mu2l/sin(psi1),  mu1j, mu2j, mubj, d1, d2, thickn_b) ;
   /* debugging 
   printf("returnwaarde mdown is %g \n", tmpnep);  /* */

         result *= tmpnep ;

   /* debugging 
   printf("result van sij_lambda_down is %g \n", result);    /* */
    

/************************************
*   debugging  
   printf("========================\n");
   printf("Cj is %f \n", Cj);
   printf("Jjlambda is %e \n", Jjlambda);
   printf("omegaj is %e \n",omegaj);
   printf("gj is %e \n",gj);
   printf("Jilambda is %e \n", Jilambda);
   printf("Jilambdaj is %e \n",Jilambdaj);
   printf("taujlambda is %e \n", taujlambda);
   printf("tauilambda is %e \n", tauilambda);
   printf("tauij os %e \n",tauij);
   printf("d1 is %g \n", d1);
   printf("d2 is %g \n", d2);
   printf("thickn_b is %g \n", thickn_b);
   printf("zi is %d \n", zi);
   printf("zj is %d \n", zj);
   printf("lambdaj is %g \n", lambdaj);
   printf("lambdai is %g \n", lambdai);
   printf("mubl %g mubj %g \n", mubl, mubj);
   printf("mu2l %g mu2j %g mu1j %g mu1i %g mu1lambda %g",
           mu2l, mu2j, mu1j, mu1i, mu1lambda);
   printf("====================\n");
    /*  */ 
return(result);
}

