/*
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*	bldatomd.c
*	builds data structure
* 	with data of all elements
*	i.e. Z, symbol, atomic weight
*	densitity
*
*	M.Bos
*	august 1997
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"


char reverse_lookup[32767];

void setatomtable(char *atomtable,char *name,int atomnumber) {
  reverse_lookup[(((int)name[0])<<8)+(int)name[1]]=atomnumber+1;
}

void bld_atlist(FILE *fptest, elem_list *lijst)
   {
    int i, tot_elems ;
    char regel[LINELENGTH];
    char blankline[LINELENGTH];
    char el_symbol[3];
    atom_dat *pntr_adat;
    int readlength, doorgaan, volgnr ;
    double atoom_gewicht, rho ;
    char  *position ;

    doorgaan=1;
    fgets(regel, LINELENGTH, fptest);
    sscanf(regel,"%d", &tot_elems);
    /* debugging
    printf("Number of elements in table is %d \n",tot_elems); /* */
    (*lijst).nr_elements=tot_elems;   
    pntr_adat= (atom_dat *) malloc(tot_elems*sizeof(atom_dat));
    (*lijst).patom_dat=pntr_adat;
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
        printf("Error in atomic data table\n");
        printf("Cannot find Z \n");
        exit(1);
     }
     position=strchr(regel, '=');
     if (position==NULL)
     {
        printf("Error in atomic data table\n");
        printf("Cannot find = \n");
        exit(1);
     }
     sscanf(position+1,"%d",&volgnr);
     /* debugging
     printf("volgnr is %d \n",volgnr); /* */
     if (volgnr !=(i+1))
     {
       printf("Error in atomic data table\n");
       printf("wrong at Z= %d \n",i+1);
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

/* put the atomic number in the structure */
    
   pntr_adat[i].z=i+1;
  
/* now get the symbol */
   
    fgets(regel,LINELENGTH,fptest);
    sscanf(regel,"%s",el_symbol);
    strcpy(pntr_adat[i].sym,el_symbol);
    /* debugging
    printf("Symbol is %s \n",el_symbol);  /* */

/* now get the atomic weight */
    fgets(regel,LINELENGTH, fptest);
    sscanf(regel,"%lf",&atoom_gewicht);
    pntr_adat[i].at_weight=atoom_gewicht;
    /* debugging
    printf("atomic weight is % g \n",atoom_gewicht);  /* */


/* now the density */
    fgets(regel,LINELENGTH, fptest);
    sscanf(regel,"%lf",&rho);
    pntr_adat[i].density=rho;
    /* debugging
    printf("density is %f \n",rho); /* */
     fgets(regel,LINELENGTH,fptest);
     readlength=strlen(regel);

     memset(blankline,32,readlength);
     if (strncmp(blankline,regel,readlength-1)!=0)
      {
        printf("No blank line found where it should be \n");
        exit(1);
      }

/* Now store element in the reverse lookup table. */
   setatomtable(lijst->inverse_lookup,pntr_adat[i].sym,i);

   }
}

void chk_atlist(elem_list *lijst)
{
    int i, z, nr_elems ;
    atom_dat *pntr_atoom ;
    char symbol[3];
    double rho, atw ;

    nr_elems= (*lijst).nr_elements;
    printf("number of elements in table is %d \n",nr_elems);
    pntr_atoom= (*lijst).patom_dat;
    for (i=0; i < nr_elems ; i++)
    {
       z=pntr_atoom[i].z;
       strcpy(symbol,pntr_atoom[i].sym);
       rho=pntr_atoom[i].density; 
       atw=pntr_atoom[i].at_weight;
       printf("Z %d  , %s : atw %f , rho %f \n",z, symbol, atw, rho);
    }
}

char *get_asymb(int z, elem_list *lijst)
{
   atom_dat *pntr_at;

   pntr_at= (*lijst).patom_dat;
   return(pntr_at[z-1].sym);
}

/************************************************************
uitgezet ivm gebruik hashtable

int get_z(char *symbol, elem_list *lijst)
{
   int i;
   atom_dat *pntr_at ;

   pntr_at= (*lijst).patom_dat;
   i=0;
   while ( strcmp(symbol,pntr_at[i].sym) !=0) i++ ;
   return(i+1);
} 
*/



/* int get_z(char *symbol, elem_list *lijst)
{
	return reverse_lookup[(((int)symbol[0])<<8)+(int)symbol[1]];
}
*/
double get_rho(int z, elem_list *lijst)
{
   atom_dat *pntr_at;
   pntr_at= (*lijst).patom_dat;
   return(pntr_at[z-1].density);
}

double get_atw(int z, elem_list *lijst)
{
   atom_dat *pntr_at;
   pntr_at= (*lijst).patom_dat;
   return(pntr_at[z-1].at_weight);
}
      
