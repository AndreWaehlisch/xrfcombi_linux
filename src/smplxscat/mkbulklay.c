/*******************************************************
*
*  COPYRIGHT (C) M.BOS 1999, FOR DETAILS SEE COPYING
*
*
*  mkbulklay.c
*
*  setups structure for a layer consisting of a single
*  element with infinite thickess
*  routine returns a pointer to this structure
*
*  M.Bos
*  oct 1997
*
*********************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

layer *mkbulklay(char *elsym)
{
   layer *pbulk ;
   elem *pelem;

   /*  create storage for layerdata and element dat */
 
   pbulk= (layer *)malloc(sizeof(layer));
   pelem=(elem *)malloc(sizeof(elem));
   strcpy((*pelem).name,elsym) ;
   (*pelem).wfract = 1.0;
    (*pbulk).nr_elements=1 ;
    (*pbulk).massthickness = 100.0 ; /* lekker dik ! */
    (*pbulk).elementen=pelem;
 
    return(pbulk);
}

  

