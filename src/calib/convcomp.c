/******************************************
*
*   COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*   convcomp.c
*
*   converts composition data in terms of compounds
*   into composition data in terms of elements
*
*   M.Bos
*   okt 1997
*
*****************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

extern FILE *own_fopen(char *);
extern elem_list lijst_els;

void convcomp(multicomplayer *multicomplay, multilayer *multil)
{
  int i,j,k, j_elem;
  int nr_layers, nr_compounds, nr_elements, z_el, nr_atoms;
  double  molwt, atfract[100], wfract, atwt;
  complayer *compslay;
  elem *reeks;
  elem_cnt *pt_elem;
  char *comp_name ;
  layer **setlagen;

  nr_layers=get_nr_clay(multicomplay);
  (*multil).nr_layers=nr_layers;
  setlagen=(layer **)malloc(nr_layers*sizeof(layer));
  if (setlagen==NULL)
   {
         printf("Cannot malloc mem for setlagen in convcomp \n");
         exit(1);
   }
  for (i=0; i < nr_layers; i++)
  {
    for (j=0;j <100; j++)
      atfract[j]=0.0;
    compslay=((*multicomplay).lagen[i]);
    nr_compounds=(*compslay).nr_compounds;
    /* debugging
    printf("laag %d nr_compounds is %d\n",i, nr_compounds); /* */
    for (j=0;j < nr_compounds;j++)
    {
      comp_name=(*compslay).comps[j].name;
      /* debugging
      printf("Compound %s \n",comp_name); /* */
      wfract=(*compslay).comps[j].wfract;
      /* debugging
      printf("weightfraction %g \n",wfract); /* */
      nr_elements=(*compslay).comps[j].nr_elems;
      /* debugging
      printf("nr_elements %d \n",nr_elements); /* */
      pt_elem=((*compslay).comps[j]).elementen;
      /* debugging
      printf("na assignment pt_elem\n"); /* */
      molwt=0.0;
      for (k=0; k < nr_elements; k++)
      {
        /* debugging
        printf("el name is %s \n",pt_elem[k].name);   /* debugging */
        z_el=get_z(pt_elem[k].name, &lijst_els);
        /* debugging
        printf("Z is %d \n",z_el); /* */
        nr_atoms=pt_elem[k].count;
        /* debugging
        printf("nr_atoms is %d \n",nr_atoms); /* */
        atwt=get_atw(z_el, &lijst_els);
        /* debugging
        printf("atwt is %g \n",atwt); /* */
        molwt += atwt*nr_atoms;
      }
      /* debugging
      printf("Na loop calcn molwt is %g \n",molwt); /* */
      for (k=0; k < nr_elements; k++)
      {
        z_el=get_z(pt_elem[k].name, &lijst_els);
        nr_atoms=pt_elem[k].count;
        atwt=get_atw(z_el, &lijst_els);
        atfract[z_el] += nr_atoms*wfract*atwt/molwt;
      }
    }
    /* now count the # of elements */
    j_elem=0;
    for (k=0; k < 100; k++)
       if (atfract[k] != 0.0) j_elem++;
  
    /* create the memory space for the elemental data */

    reeks=(elem *) malloc(j_elem*sizeof(elem));
    /* debugging
    printf("Na alloc reeks \n"); /* */
    if (reeks==NULL)
    {
      printf("Cannot malloc mem for reeks in convcomp \n");
      exit(1);
    }

    /* now start to fill in data structure with elem. data */

   j_elem=0;
   for (k=1; k < 100; k++)
   {
      if (atfract[k] != 0.0)
      {
        /* printf("in all elements loop something found for k is %d\n",k); */
         reeks[j_elem].wfract =atfract[k];
         strcpy(reeks[j_elem].name, get_asymb(k, &lijst_els));
         j_elem++;
      }
    }
   setlagen[i] = (layer *) malloc(sizeof(layer));
   setlagen[i]->nr_elements= j_elem;
   setlagen[i]->elementen = reeks;
   setlagen[i]->massthickness=(*compslay).massthickness;
  }
  (*multil).lagen=setlagen;
/* debugging
checkdata(multil); /* */
}

