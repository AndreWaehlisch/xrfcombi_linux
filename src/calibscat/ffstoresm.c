/*
*    file ffstore.c
*    routine to fill datastructure
*    with table data formfactors
*
*    M.Bos
*    Dec 6th 2003
*
*
*    Dec 14th 2003
*    changed code to accomodate
*    sm_tables with slightly different
*    TITLE = heading
*    
*    Dec 27th 2003
*    changed to store formfactor instead of cross section
*    and changed calcff2 to use formfactor
*     
*/
            
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"


/* globals */
 
extern double cross [100][60][100] ;
extern double energies [100][60] ;
extern double thetas [100][100];


int storeffs(void)
{
	int z ;
	FILE *fp ;
	char full_line[LINELENGTH];
	char filename[80];
	char *chknull;
	char partname[]="_cs0sl_mf" ;
	int  found, fndend, fndblock ;
	char *pos1, *pos2;
    char *lineptr;
	double curr_energy;
	int energy_cnt, angle_cnt;
	double angle_lo, angle_curr, scatter;
	double  dum, formfactor ;

	/* loop over 99 elements in the datatables */
	for (z=1; z < 100; z++)
	{
		/* construct filename */
		sprintf(filename, "./scattercoefs/%03d%s", z, partname);
		if ((fp=fopen(filename,"r"))==NULL)
		{
			printf("Cannot open %s\n",filename);
			exit(1);
		}
		energy_cnt=1;

lus:
        fndblock=0;
	

    /* search for BLOCK */
		do
		{
			chknull = fgets(full_line, LINELENGTH, fp);
			if (chknull==NULL)
			{ 
				fndend=1;
				/* printf("end of file?? \n"); */
				goto finale;
			}
			/* debugging 
			puts(full_line) ; /* */
		}
		while ((pos1=strstr(full_line,"*BLOCK:"))==NULL) ;
      
		lineptr=full_line ;
		/* strip string to find energy */


        /*   title line in _sm files */
		lineptr=strtok(lineptr,":");  /*  BLOCK*/
		pos2=strtok(NULL, "k" );  /* 0.05430 energy in keV */
		/* debugging 
		printf("%s \n", pos2) ; /* */
        
		curr_energy=strtod(pos2,NULL);
        /* debugging 
		printf("Energy %g \n", curr_energy);  /*  */

        energies[z][energy_cnt]= curr_energy;

		/* now find line starting with THETA */
		do
		{
			chknull=fgets(full_line, LINELENGTH, fp);
			if (chknull==NULL)
			{
				printf("Cannot find THETA line \n");
			    exit(1);
			}
		} while ((pos1=strstr(full_line,"THETA"))==NULL);

	    /* found THETA line, now get angle data up to 180 degrees */
		found=0;
		angle_lo=0.0;
		angle_curr=0.0;
		angle_cnt=1;
		do
		{
			chknull=fgets(full_line, LINELENGTH, fp);
			if (chknull==NULL)
			{
				printf("Error in reading line with angle data\n");
				exit(1);
			}
			sscanf(full_line,"%lf %lf %lf %lf",&angle_curr, &scatter, &dum, &formfactor);
			if (angle_curr > 179.999)
			{
				found=1;
			}
            thetas[z][angle_cnt]=angle_curr;
			cross[z][energy_cnt][angle_cnt++]=formfactor;
		} while (found==0);
		energy_cnt++;
		goto lus ;
finale:
	fclose(fp);
	energy_cnt=1;
	}
	return (0);
}



double calcff2(int z, double theta, double energy)
{
	double energy_lo, energy_hi, curr_energy ;
	double angle_lo, angle_hi, curr_angle;
	double scatterfactor, scatter_lo, scatter_hi;
	double scatter_mean_lo, scatter_mean_hi;
	int i,j, found, error, nr_energy, nr_theta ;
    double costheta ;

     /* find energy bounds */

     found=0;
	 error=0;
	 i=1; 
	 do
	 {
		 curr_energy=energies[z][i];
		 if (curr_energy > energy)
		 {
			 energy_hi=curr_energy;
			 found=1;
		 }
		 else 
		 {
			 energy_lo=curr_energy;
			 i++;
			 if (i > NR_ENERGIES)
			 {
				 error=1;
				 found=1;
			 }
		 }
	 } while (found==0);
	 nr_energy=i;

	 if (error)
	 {
		 printf("Energy bounds not found \n");
		 printf("Last energy in table is %g \n", energies[z][i]);
		 exit(1);
	 }

	 /* find theta bounds */
	 found=0;
	 error=0;
	 j=1;
	 do
	 {
		 curr_angle=thetas[z][j];
		 if (curr_angle >= theta)   /* PA0MBO 22 dec 2003  */
		 {
			 angle_hi=curr_angle;
			 found=1;
		 }
		 else
		 {
			 angle_lo=curr_angle;
			 j++;
			 if (j > NR_THETAS)
			 {
				 error=1;
				 found=1;
			 }
		 }
	 } while (found==0);
	 nr_theta=j;

	 if (error)
	 {
		 printf("Theta bounds not found \n");
		 printf("Last theta in tables is %g\n", thetas[z][j]);
		 exit(1);
	 }

	 /* now get cross sections */
	 scatter_lo = cross[z] [nr_energy-1][nr_theta-1];
	 scatter_hi = cross[z] [nr_energy-1][nr_theta];
   /*  printf("scatter_lo %g scatter_hi %g nr_theta %d \n", scatter_lo, scatter_hi, nr_theta);  /*  */
     if (nr_theta==1)
		 scatter_mean_lo=scatter_hi ;
	 else
	   scatter_mean_lo = scatter_lo +(theta-angle_lo)/(angle_hi-angle_lo)*(scatter_hi-scatter_lo);
/*  printf("angle_hi %g angle_lo %g scatter_mean_lo %g\n", angle_hi, angle_lo, scatter_mean_lo);  /*  */
     scatter_lo = cross[z] [nr_energy][nr_theta-1];
	 scatter_hi = cross[z] [nr_energy][nr_theta];
  /* printf("scatter_lo %g scatter_hi %g \n", scatter_lo, scatter_hi);  /*  */
	 if (nr_theta==1)
		 scatter_mean_hi = scatter_hi ;
     else
        scatter_mean_hi = scatter_lo +(theta-angle_lo)/(angle_hi-angle_lo)*(scatter_hi-scatter_lo); 
	 scatterfactor=scatter_mean_lo+(energy-energy_lo)/(energy_hi-energy_lo)*
		           (scatter_mean_hi-scatter_mean_lo);
	 /* debugging 
	 printf("energy_lo %g energy_hi %g\n", energy_lo, energy_hi);  
	 printf("ff uit interpolatie tabel %g \n", scatterfactor); /* */
	 scatterfactor *= scatterfactor ;
	 scatterfactor *= R0*R0/2.0;
         costheta =cos( (theta/180)*PI);
  	 scatterfactor *= (1.0+ costheta*costheta);
     /* debugging 
	 printf("cross section from formfactor is  (cm2/atom ) %g theta %g\n", 
		   scatterfactor, theta); /* */

	 return(scatterfactor);
}
