/*
*   file monosim.c
*
*   calculates contribution
*   of compton scattering followed 
*   by  photoelectric event
*  
*   based on Tirao's paper
*   
*
*   double integral over phi and cosa
*   scattering angle theta calculated from phi and alfa
*   labdastar (belonging to scattered photon)
*   calculated from  theta and compton equation
*
*   M.Bos
*
*   Jan 17th 2004
*
*
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

int main(void)
{
  multilayer sample ;
  multilayer bulksample ;
  layer *ptop, *pbulk ;
  double primint, secint, totint, comptonint, primint_bulk, comptonint_bulk, totint_bulk;
  double  totcomptonint, totcomptonint_bulk, totprimint, totprimint_bulk ;
  double totsecint ;
  double raylfotint, raylfotint_bulk, totraylfotint, totraylfotint_bulk;
  double fotraylint, fotraylint_bulk, totfotraylint, totfotraylint_bulk;

  struct spectrum_pair *pntrspectrum ;
  double labda_in, excitation_energy ;
  int i,j ;
  char el[3];
  char lin[6];
  int z_filt;
  double mth_filt, kev ;
  FILE *fpresult ;
 
  fpresult=fopen("uitvoer","w");
  init_tables();
  set_spectrom();
  storeffs();
  def_sample(&sample);
  ptop=get_layer(&sample, 0);
 
  printf("nr lines is %d\n", nr_of_meas_lines);

  
  for (j=0; j < nr_of_meas_lines; j++)
  {
    f2cnt=0;   /* pa0mbo tbv registratie functiecalls f2 in integratie */
    strcpy(el, get_asymb((lines_to_meas[j]).z, &lijst_els));
    strcpy(lin, (lines_to_meas[j]).line);
    z_filt=(lines_to_meas[j]).filt_z;
    mth_filt=(lines_to_meas[j]).mth_filt;
    kev=(lines_to_meas[j]).kev;
    pntrspectrum=gen_spec(anode_el, take_off, berryliumd, kev, 
                          PNTS_SPECTRUM, z_filt, mth_filt);
    pbulk=mkbulklay(el);
    bulksample.nr_layers=1;
    bulksample.lagen=&pbulk;

    totint=0.0;
    totsecint=0.0;
    totint_bulk=0.0;
    totcomptonint=0.0;
    totcomptonint_bulk=0.0;
    comptonint=0.0 ;
    comptonint_bulk=0.0 ;
    totprimint=0.0;
    totprimint_bulk=0.0;
    totraylfotint =0.0;
    totraylfotint_bulk=0.0;
    raylfotint=0.0;
    raylfotint_bulk =0.0;
    fotraylint =0.0;
    fotraylint_bulk = 0.0;
    totfotraylint = 0.0;
    totfotraylint_bulk = 0.0;
   


  	  printf("Give excitation energy keV\n");
           scanf("%lf", &excitation_energy);
	  labda_in = 12.396/excitation_energy ;
	  i=0; 
	  (pntrspectrum[i]).intens=1;  /* */


      primint= pilabda(el,lin,labda_in, ptop)*(pntrspectrum[i]).intens;
      primint_bulk = pilabda(el,lin,labda_in, pbulk)*(pntrspectrum[i]).intens;
      secint=intrasecflu(el,lin,labda_in, &sample, 0, primint) ; 
      comptonint=comptisesr(el,lin, labda_in, ptop,i, 0) *(pntrspectrum[i]).intens; 
      comptonint_bulk=comptisesr(el, lin,labda_in, pbulk,i, 1)*(pntrspectrum[i]).intens; 
      raylfotint=raylisesr(el,lin, labda_in, ptop,i, 0) *(pntrspectrum[i]).intens; 
      raylfotint_bulk=raylisesr(el, lin,labda_in, pbulk,i, 1)*(pntrspectrum[i]).intens; 
      fotraylint=ispfr(el,lin, labda_in, ptop,i, 0) *(pntrspectrum[i]).intens; 
      fotraylint_bulk=ispfr(el, lin,labda_in, pbulk,i, 1)*(pntrspectrum[i]).intens; 
      /* debugging  
      printf("pnt %d %g %g %g %g %g\n", i,primint, primint_bulk, secint,
              comptonint,  comptonint_bulk); /* */

      totprimint += primint;
      totsecint += secint;
      totcomptonint +=comptonint;
      totraylfotint += raylfotint;
      totfotraylint += fotraylint;
      totint += primint+secint+comptonint+raylfotint+ fotraylint;

      totprimint_bulk += primint_bulk; 
      totcomptonint_bulk+=comptonint_bulk;
      totraylfotint_bulk+=raylfotint_bulk;
      totfotraylint_bulk+=fotraylint_bulk;
      totint_bulk +=primint_bulk+comptonint_bulk+raylfotint_bulk+fotraylint_bulk;


    printf("%s %s totint %g totint_bulk %g totprimint %g totprimintbulk %g"
                     "totsecint %g totcomptonint %g totcomptonintbulk %g\n"
               "totraylfotint %g totraylfotint_bulk %g \n"
               "totfotraylint %g totfotraylint_bulk %g  \n",
            el, lin, totint, totint_bulk, totprimint, totprimint_bulk,
                         totsecint, totcomptonint, totcomptonint_bulk,
                  totraylfotint, totraylfotint_bulk, totfotraylint, totfotraylint_bulk);
	printf("rxi old is %g rxi scatter is %g \n", 
		(totprimint+totsecint)/(totint_bulk-totcomptonint_bulk-totraylfotint_bulk-
		   totfotraylint_bulk), totint/totint_bulk); 
/*	printf("f2 functiecalls cnt is %d\n", f2cnt);  */
    fprintf(fpresult,"%s %s %e\n",el,lin,totint/totint_bulk);
 }  
  fclose(fpresult);
 return 0;
}
         




