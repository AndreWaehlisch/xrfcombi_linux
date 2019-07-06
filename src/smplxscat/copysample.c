/***************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*    copysample.c
*
*    routine to build a duplicate
*    sample structure  
*
*    M.Bos
*    februari 1998
*
***********************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

void copysample(multilayer *source, multilayer *dest)
{
  int i,j ;
  int nr_layers, nr_elements ;
  layer **setlagen;
  elem *reeks ;

  nr_layers= (*source).nr_layers;
  (*dest).nr_layers=nr_layers ;
  setlagen=(layer **)malloc(nr_layers*sizeof(layer));
  (*dest).lagen=setlagen ;
  for (i=0; i < nr_layers ; i++)
  {
   (*dest).lagen[i]= (layer *)malloc(sizeof(layer));
   nr_elements = (*(*source).lagen[i]).nr_elements;
   (*(*dest).lagen[i]).nr_elements=nr_elements;
   (*(*dest).lagen[i]).massthickness=(*(*source).lagen[i]).massthickness;
   reeks = (elem *) malloc( nr_elements*sizeof(elem));
   (*(*dest).lagen[i]).elementen=reeks;
   for (j=0; j < nr_elements ; j++)
   {   
     strcpy( ((*(*dest).lagen[i]).elementen[j]).name, 
             ((*(*source).lagen[i]).elementen[j]).name);
     ((*(*dest).lagen[i]).elementen[j]).wfract=
         ((*(*source).lagen[i]).elementen[j]).wfract;
   }
  }
  /* debugging 
  checkdata(source);  /* */
}
             

