/*
*
*    file calcsxz.c
*
*    calculation of incoherent atomic scattering
*    function based on thomas-fermi model
*
*    eqns from J.E. Fernandez, V.G. Molinari and M. Sumini
*    Adv. X-ray Anal. 33 (1990) 553
*
*    coded by M.Bos
*    Jan 7th 2004
*
*/
#include <stdio.h>
#include <math.h>
#include "xrfluor.h"

#define E_electron 510.98                 /* energy electron at rest keV */

double calcsxz(int z, double theta, double energy)
{
	double labda, costheta, vgreek, result;
	double kkn, sigma, labda_prime, dlabda, term2, sfactor;
    double Ec, eps, gamma, gamma_prime ;
	labda = 12.396/energy ;
/*	printf("labda is %g \n", labda) ;/* */
	costheta = cos (theta*PI/180.0);
/*	printf("costheta %g \n", costheta) ; /* */

	/* system fernandez */ 
	vgreek = (2.0/3.0)*(137.0/pow((double) z,(2.0/3.0)));	
/*	printf("vgreek1 %g \n", vgreek); /* */
	vgreek *= (Compton/labda);
/*	printf("vgreek2 %g \n", vgreek); /* */
	vgreek *= sqrt( (1.0-costheta)/2.0);
/*	printf("vgreek is %g \n",vgreek); /* */
    sfactor = 1.0 - exp( -4.88* pow(vgreek, 0.856));
/*	printf("Sfactor in calcsxz %g theta %g labda %g\n", sfactor, theta, labda); /* */
	sigma = R0*R0/2 ;
	dlabda = (1.0-costheta)*Compton ;
/*	printf("in calcsxz dlabda is %g \n", dlabda);  /*  */
	labda_prime = labda + dlabda ;

	/* system fernandez 
	term2 = labda_prime/labda ;
	kkn= term2*term2*(1/term2 + term2  + (dlabda/Compton)*(dlabda/Compton -2.0));
    printf("kkn fernandez %g  kkn*sigma %g \n", kkn, kkn*sigma);  /* */


	/* system Janos Kirz */
	term2=  (1+(energy/E_electron)*(1-costheta));
	kkn = (1+costheta*costheta)/(term2*term2);
/*	printf("kkn Janos Kirz is %g kkn*sigma %g \n", kkn, kkn*sigma); /*  */

	/* system derek paul web 
    gamma_prime = (12.396/labda_prime)/E_electron ;
	gamma= energy/E_electron ;
    term2= gamma_prime/gamma ;
	kkn= term2*term2*(1/term2 + term2 - sin(theta*PI/180)*sin(theta*PI/180));
    printf("kkn derek paul %g kkn*sigma %g \n", kkn, kkn*sigma);

	/*  sampling final state web 

	Ec=  energy* (E_electron/(E_electron+energy*(1-costheta)));
    eps=Ec/energy;
	kkn= (E_electron/energy)*(1/eps + eps)*(1- (eps*sin(theta*PI/180)*sin(theta*PI/180))/(1+eps*eps));
	printf("kkn sampling  is %g kkn*sigma %g \n", kkn, kkn*sigma); /*  */
	result = sigma*kkn*sfactor ;
/*	printf("labda %g labdaprime %g result %g\n", labda, labda_prime, result); /* */
	return result;
}

  