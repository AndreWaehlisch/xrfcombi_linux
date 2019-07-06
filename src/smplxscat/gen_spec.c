/******************************************************
*
*    COPYRIGHT (C) M.BOS 1999 , for details see copying
*
*
*    gen_spec.c
*
*    generates X-ray tube spectrum based
*    on Pella's algorithms
*    P.A.Pella, L. Feng, J.A.Small
*    X-ray spectrom. 14 (1985) 125
*
*    coded by
*    M.Bos
*    oct 1997
*
*    adjusted characteristic line treatment
*    Oct 1998
*
*    attempt to adjust influence charac. lines from tube
*    spectrum to correct lambda-range, i.e. not extending
*    to longer wavelengths.
*    Feb 1999
*
*    charac.line intensity spread over 2 spectrumpoints
*    jan 25th 2003
*
*
*    corrected tube filter code _ March 2004
********************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

extern elem_list lijst_els;
extern emudata abstabel ;
extern double pellas_consts[4][3];
struct spectrum_pair  *gen_spec(char *anode_el, double take_off, double wind_thickn,
              double kv, int nr_points, int z_filt, double mth_filt)
{
   int i, z_anode ;
   double xsi, psi, labd0, labda ;
   double l_intv ;                   /* lambda interval */
   double mu_target, f, c, wabs; 
   double lka, lkb1, la12, lb1 ;
   double ratioka, ratiokb1, ratiola12, ratiolb1, u0 ;
   double tmpintens1, tmpintens2, tmpintens3, tmpintens4;
   FILE *fpspectrum;
   struct spectrum_pair *spectrum ;
   int flagka, flagkb1, flagla12, flaglb1;
   int nrka, nrkb, nrla12, nrlb1 ;
   double t1, t2, t3 ;
   double  mu_filt;
   int filt_on ;

   flagka=flagkb1=flagla12=flaglb1=0;

  if ((fpspectrum=fopen("spectrum.dat","w"))==NULL)
   {
      printf("Cannot open file spectrum.dat \n");
      exit(1);
   }  /* */

   if ( mth_filt ==0.0)
      filt_on=0;
   else
        filt_on=1;
   

   z_anode=get_z(anode_el, &lijst_els);
   labd0 = 12.396/kv ;
   l_intv = (LAMBDA_HIGHEST-labd0)/nr_points;
   spectrum=  malloc(nr_points*sizeof( struct spectrum_pair));
   psi = (take_off/360.0)*2*PI ;
   lka=characwl("ka", anode_el, &lijst_els);
   if (lka < labd0) flagka=1;
   lkb1=characwl("kb1", anode_el, &lijst_els);
   if (lkb1 < labd0) flagkb1=1;
   la12=(characwl("la1", anode_el, &lijst_els)+
         characwl("la2", anode_el, &lijst_els))/2.0 ;
   if (la12 < labd0) flagla12=1;
   lb1=characwl("lb1", anode_el, &lijst_els);
   if (lb1 < labd0) flaglb1=1;

 /*  printf("in gen_spec flagka %d flagkb1 %d flagla12 %d flaglb1 %d \n",
           flagka, flagkb1, flagla12, flaglb1); /* */ 

   for (i=0; i < nr_points; i++)
   {
      if (i==0)
        labda=labd0*1.001;    /* to make sure nr 1 != 0.0 */
      else
        labda= labd0+ i*l_intv;
      
      (spectrum[i]).lambda=labda;
      mu_target=calc_muilambda(z_anode, 12.396/labda, &abstabel);
      xsi= (1.0/pow(labd0,1.65) - 1.0/pow(labda, 1.65))*(mu_target/sin(psi));
      c = 1.0 + 1.0/(1.0 +2.56E-3*z_anode*z_anode);
      c= c/((1.0+2.56E3*labd0/(z_anode*z_anode))*(0.25*xsi+1.0E4));
      f=1.0/((1.0+c*xsi)*(1.0+c*xsi));
      wabs= exp(-0.35*pow(labda,2.86)*wind_thickn);
      (spectrum[i]).intens=
            2.72E-6*z_anode*(labda/labd0-1.0)*f*wabs/(labda*labda);
    }
    if (flagka ==0)
    {
          nrka=(lka-labd0)/l_intv;
          u0=lka/labd0 ;
          t1= (u0-1.0)/(1.17*u0+3.2);
          t1=0.5*t1*t1;
          t1=exp(-t1);
          t2=pellas_consts[0][0]/(pellas_consts[0][1]+
             z_anode*z_anode*z_anode*z_anode)+pellas_consts[0][2];
          t3=u0*log(u0)/(u0-1.0)-1.0;
          ratioka=t1*t2*t3;
          tmpintens1= (spectrum[nrka]).intens*ratioka/(2*l_intv);
          (spectrum[nrka]).intens +=tmpintens1;
          (spectrum[nrka+1]).intens +=tmpintens1;
       }
      if ( flagkb1 ==0)
      {
          nrkb=(lkb1-labd0)/l_intv;
          u0=lkb1/labd0 ;
          t1= (u0-1.0)/(1.17*u0+3.2);
          t1=0.5*t1*t1;
          t1=exp(-t1);
          t2=pellas_consts[1][0]/(pellas_consts[1][1]+
             z_anode*z_anode*z_anode*z_anode)+pellas_consts[1][2];
          t3=u0*log(u0)/(u0-1.0)-1.0;
          ratiokb1=t1*t2*t3;
          tmpintens2= (spectrum[nrkb]).intens*ratiokb1/(2*l_intv);
          (spectrum[nrkb]).intens +=tmpintens2;
          (spectrum[nrkb+1]).intens +=tmpintens2;
       }
      if ( flagla12 ==0)
      {
         nrla12=(la12-labd0)/l_intv;
          u0=la12/labd0 ;
          t1= (u0-1.0)/(1.17*u0+3.2);
          t1=0.5*t1*t1;
          t1=exp(-t1);
          t2=pellas_consts[2][0]/(pellas_consts[2][1]+
             z_anode*z_anode*z_anode*z_anode)+pellas_consts[2][2];
          t3=u0*log(u0)/(u0-1.0)-1.0;
          ratiola12=t1*t2*t3;
          tmpintens3= (spectrum[nrla12]).intens*ratiola12/(2*l_intv);
          (spectrum[nrla12]).intens +=tmpintens3;
          (spectrum[nrla12+1]).intens +=tmpintens3;
       }
      if ( flaglb1 ==0)
      {
          nrlb1=(lb1-labd0)/l_intv;
          u0=lb1/labd0 ;
          t1= (u0-1.0)/(1.17*u0+3.2);
          t1=0.5*t1*t1;
          t1=exp(-t1);
          t2=pellas_consts[3][0]/(pellas_consts[3][1]+
             z_anode*z_anode*z_anode*z_anode)+pellas_consts[3][2];
          t3=u0*log(u0)/(u0-1.0)-1.0;
          ratiolb1=t1*t2*t3;
          tmpintens4= (spectrum[nrlb1]).intens*ratiolb1/(2*l_intv);
          (spectrum[nrlb1]).intens +=tmpintens4;
          (spectrum[nrlb1+1]).intens +=tmpintens4;
       }
       if (filt_on)
       {
          for (i=0; i < nr_points; i++)
          {
            labda=(spectrum[i]).lambda ;
          mu_filt=calc_muilambda(z_filt, 12.396/labda, &abstabel);
          (spectrum[i]).intens *= exp(-mu_filt*mth_filt);
          }
       } 

      for (i=0; i < nr_points; i++)
      {
          fprintf(fpspectrum,"%e  %e \n",(spectrum[i]).lambda,
              (spectrum[i]).intens); /* */

        /*   fprintf(fpspectrum,"xsi %e, c %e , f %e , wabs  %e , spec %e \n",
               xsi, c, f, wabs, spectrum[i]); */
      }
  fclose(fpspectrum); 
 /*  printf("net voor return in gen-spec \n"); /* */
   return(spectrum);
}

