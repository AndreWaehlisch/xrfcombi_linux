/*
*    file raylcrossmix.c
*
*    calculates rayleigh scattering crossection
*    for a sample with know composition
*    using scattering angle theta
*    and energy of primary radiation 
*
*    M. Bos Nov 28th 2003
*
*    input parameters
*    -  double energy : energy of radiation that is scattered
*    -  layer * laag  : pointer to layer structure of sample layer
*    -  eleme_list *lijst : pointer to structure with data of elements
*    -  double theta : scattering angle
*
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "xrfluor.h"

double raycrossmix(double energy, layer *laag, elem_list *lijst, double theta)
{
  int i, z, nr_elems;
  double wfraction, crossE, crosstmp, atw ;
  
  crossE=0.0 ;
  /* debugging 
  printf("in raycrossmix energy is %g\n", energy); /* */
  nr_elems=(*laag).nr_elements ;
  for (i=0; i < nr_elems;i++)
  {
    z=get_z(((*laag).elementen[i]).name,lijst);
	/* debugging */

	if (z> 92)
	{
		printf("unknown elemen %d nr elements was %d \n", z);
	} /* */

		
	atw=get_atw(z, lijst);
    wfraction=get_weightfract(laag,i);
    /* debugging 
    printf("in raycrossmix z is %d, theta   is %g , energy is %g\n",
      z, theta, energy); /* */ 
      crosstmp=calcff2(z,theta,energy);  

 /* debugging */
 /* printf("wt %g crosstmp %g z %d theta %g energy %g\n",  
	  wfraction, crosstmp, z, theta, energy); /* */
    crosstmp*=Na/atw; /* from barn to cm2/g /*  */
 /* printf("crosstmp %g theta %g\n",crosstmp, theta);  /*  */
    crossE += wfraction*crosstmp;
  }
 
  return(crossE);
}
 
