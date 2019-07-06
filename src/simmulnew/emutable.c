/***************************************
*
*   emutable.c
*
*   routines to build and access
*   tabel with energy/massabsorption coefs
*
*   M.Bos
*   august 1997
*
*****************************************/
#include <stdio.h>
#include  <math.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

   
void bldemutable(FILE *fptest, emudata *abstable)
{
   int i,j ;
   char regel[LINELENGTH];
   char blankline[LINELENGTH];
   int tot_elems, volgnr ;
   atom **set_elementen ;
   emupair *reeks ;
   char *position ;
   int readlength, paircount, doorgaan, verder;
   double coefread;
   double tmpstorecoefs[MAXCOEFS];
   char *startpntr, **endpntr ;
   char *eindecoef ;  

   endpntr=&eindecoef;
   doorgaan =1;
   verder=1;
   fgets(regel,LINELENGTH,fptest); 
   sscanf(regel,"%d", &tot_elems);
   /* debugging 
   printf("Number of elements in table is %d\n",tot_elems); /* */
   set_elementen= (atom **) malloc(tot_elems*sizeof(atom));
     fgets(regel,LINELENGTH,fptest);
     readlength=strlen(regel);
     /* debugging 
     printf("readlength = %d\n",readlength);  /* */
     memset(blankline,32,readlength);
     blankline[readlength-1]='\0' ;
     if (strncmp(blankline,regel,readlength-1)!=0)
      {
        printf("No blank line found where it should be \n");
        printf("Offending line is %s \n",regel);
        printf("blank line is %s\n",blankline);
        exit(1);
      }
   for (i=0; i < tot_elems; i++)
   {
     fgets(regel,LINELENGTH,fptest); /* line with Z = .. */
     position=strchr(regel, 'Z');
     if (position==NULL)
     {
        printf("Error in mass-absorption data table\n");
        printf("Cannot find Z \n");
        exit(1);
     }
     position=strchr(regel, '=');
     if (position==NULL)
     {
        printf("Error in mass-absorption data table\n");
        printf("Cannot find = \n");
        exit(1);
     }
     sscanf(position+1,"%d",&volgnr);
     /* debugging 
     printf("volgnr is %d \n",volgnr); /* */
     if (volgnr !=(i+1))
     {
       printf("Error in mass-absorption data table\n");
       printf("Cannot wrong at Z= %d \n",i+1);
       exit(1);
     }
     fgets(regel,LINELENGTH,fptest);
     readlength=strlen(regel);
     memset(blankline,32,readlength);
     if (strncmp(blankline,regel,readlength-1)!=0)
      {
        printf("No blank line found where it should be \n");
        exit(1);
      }

    /* now get coeffs */
    
    paircount=0;
    do
    {
     fgets(regel,LINELENGTH,fptest);
     readlength=strlen(regel);
     memset(blankline,32,readlength);
     if (strncmp(blankline,regel,readlength-1)!=0)  /* line with coeffs */
      {
          doorgaan=1;
          startpntr=regel;
          do
          {
               coefread=strtod(startpntr, endpntr);
               tmpstorecoefs[paircount]=coefread;
               /* debugging 
               printf(" coef is %f \n",coefread); /* */
               paircount++;
               startpntr=eindecoef;
          } while ( *endpntr < regel+readlength-4 ); /* max 4 spaces at eol */
        if ( (paircount % 2)!= 0)
         {
 printf("uneven number of coefs \n");
            exit(1);
         }
      }
     else
         doorgaan=0;
     } while (doorgaan);
   reeks= (emupair *)malloc((paircount/2)*sizeof(emupair));
   for (j=0; j < paircount/2 ; j++)
   {
        reeks[j].energy = tmpstorecoefs[j*2];
        reeks[j].abscoef = tmpstorecoefs[j*2 + 1];
   }
   set_elementen[i]= (atom *) malloc( sizeof(atom));
   set_elementen[i]->nr_emupairs=paircount/2;
   set_elementen[i]->paren=reeks;
   set_elementen[i]->z=volgnr;
  }
  (*abstable).at_nrs=tot_elems;
  (*abstable).atomcoefs=set_elementen;
}


void chkdata( emudata *pemudata)
{
  int aant_atnrs, i, j ;
  atom *patom;
  emupair *pemupair;
  int nr_pairs, z;

  aant_atnrs = (*pemudata).at_nrs;
  printf("Totaal aantal entries in tabel is %d \n",aant_atnrs);
  for (i=0; i < aant_atnrs; i++)
  {
    patom=(*pemudata).atomcoefs[i];
    z= (*patom).z ;
    printf("Z is %d \n",z);
    nr_pairs= (*patom).nr_emupairs;
    printf("Aantal paren voor Z= %d is %d\n",z, nr_pairs);
    pemupair= (*patom).paren;
    for (j=0; j < nr_pairs; j++)
    {
       printf(" E= %f  mu is  %f \n",pemupair[j].energy, pemupair[j].abscoef);
    }

  }
}

double calc_muilambda(int z, double kev, emudata *pemudata)
{
  int i ;
  double kevbig, kevsml, mu, musml, mubig ;
  double A, slope ;
  int nr_pairs  ;
  emupair *pemus;
  atom *patom;
  /* first check z */
  if ( ( z < 1) || ( z> HIGHEST_Z))
  {
       printf("Atomic number %d is out of range \n",z);
       exit(1);
  } 
  patom= (*pemudata).atomcoefs[z-1];
  nr_pairs=(*patom).nr_emupairs;
  pemus= (*patom).paren;
  i=0;
  while ( pemus[i].energy < kev )
  {
     i++;
  }
  mubig=pemus[i].abscoef;
  kevbig=pemus[i].energy;
  musml=pemus[i-1].abscoef;
  kevsml=pemus[i-1].energy;
  slope= ( log(mubig) - log(musml))/(log(kevsml)-log(kevbig));
  A= log(mubig) + slope*log(kevbig);
  mu = exp(A-slope*log(kev));
  return(mu);
}


double calc_mumix(double energy, layer *laag, elem_list *lijst,
                   emudata *pemudata)
{
     int i, z, nr_elems ;
     double wfraction, mulabda, mutemp ;

     mulabda=0.0;
     nr_elems=(*laag).nr_elements ;
     for (i=0; i < nr_elems ; i++)
     {
        z=get_z( ((*laag).elementen[i]).name, lijst);
        mutemp=calc_muilambda(z, energy, pemudata);
/*        wfraction= ((*laag).elementen[i].wfract); */
       wfraction=get_weightfract(laag,i); 
/*       printf(" mutemp %e wfraction %f \n",mutemp,wfraction); */
        mulabda += mutemp*wfraction;
/*        mulabda += get_weightfract( laag, i)*calc_muilambda(z, energy, pemudata); */
      }
     return(mulabda);
}
    
