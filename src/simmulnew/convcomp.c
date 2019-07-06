/******************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
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
extern double sumfracts[MAX_LAYERS] ;

void convcomp(multicomplayer *multicomplay, multilayer *multil)
{
  int i,j,k, j_elem;
  int nr_layers, nr_compounds, nr_elements, z_el, nr_atoms;
  double molwt, atfract[100], wfract, atwt;
  complayer *compslay;
  elem *reeks;
  elem_cnt *pt_elem;
  char *comp_name ;
  layer **setlagen;

  nr_layers=get_nr_clay(multicomplay);
  (*multil).nr_layers=nr_layers;
  setlagen=(layer **)malloc(nr_layers*sizeof(layer));
  for (i=0; i < nr_layers; i++)
  {
    sumfracts[i] = 0.0;
    for (j=0;j <100; j++)
      atfract[j]=0.0;
    compslay=((*multicomplay).lagen[i]);
    nr_compounds=(*compslay).nr_compounds;
/*     printf("laag %d nr_compounds is %d\n",i, nr_compounds);  /* */
    for (j=0;j < nr_compounds;j++)
    {
      comp_name=(*compslay).comps[j].name;
/*    printf("Compound %s \n",comp_name);  /* */
      wfract=(*compslay).comps[j].wfract;
      if (wfract==0.0) wfract=1.0E-32;  /* to prevent missing of el */
      sumfracts[i] +=wfract;
/*  printf("weightfraction %g \n",wfract);  /* */
      nr_elements=(*compslay).comps[j].nr_elems;
 /*     printf("nr_elements %d \n",nr_elements); /* */
      pt_elem=((*compslay).comps[j]).elementen;
      /* printf("na assignment pt_elem\n"); */
      molwt=0.0;
      for (k=0; k < nr_elements; k++)
      {
        z_el=get_z(pt_elem[k].name, &lijst_els);
/*        printf("Z is %d \n",z_el);/* */
        nr_atoms=pt_elem[k].count;
/*       printf("nr_atoms is %d \n",nr_atoms); /* */
        atwt=get_atw(z_el, &lijst_els);
/*        printf("atwt is %g \n",atwt); /* */
        molwt += atwt*nr_atoms;
      }
/*       printf("Na loop calcn molwt is %g \n",molwt); /* */
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
/*    printf("Count of nr of elements in whole layer is %d \n", j_elem); /* */  
    /* create the memory space for the elemental data */

    reeks=(elem *) malloc(j_elem*sizeof(elem));
    /* printf("Na alloc reeks \n"); */

    /* now start to fill in data structure with elem. data */

   j_elem=0;
   for (k=1; k < 100; k++)
   {
      if (atfract[k] != 0.0)
      {
        /* printf("in all elements loop something found for k is %d\n",k); /* */
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
/*  checkdata(multil);  /* debugging */
}


void free_multil(multilayer *multil)
{
   int nr_layers, i ;
   layer **setlagen ;
 
   nr_layers= (*multil).nr_layers;
   setlagen= (*multil).lagen;
   for (i=0; i < nr_layers; i++)
   { 
      free(setlagen[i]->elementen);
      free(setlagen[i]);
   }
   free(setlagen);
}
