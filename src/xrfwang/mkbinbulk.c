/*******************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*  mkbinbulk.c
*
*  setups structure for a multilayer consisting of 
   one layer with two 
*  elements with infinite thickess
*  routine returns a pointer to this structure
*
*  M.Bos
*  jan  1998
*
*********************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

multilayer *mkbinbulk(char *elsym1, char *elsym2, double wfract1)
{
   layer *pbulk ;
   elem *pelem ;
   multilayer *pmlayer;

   /*  create storage for layerdata and element dat */
 
   pbulk= (layer *)malloc(sizeof(layer));
   pelem=(elem *)malloc(2*sizeof(elem));
   pmlayer=(multilayer *)malloc(sizeof(multilayer));
   (*pmlayer).lagen=(layer **)malloc(sizeof(layer)); 
   strcpy((pelem[0]).name,elsym1) ;
   strcpy((pelem[1]).name,elsym2) ;
   (pelem[0]).wfract = wfract1;
   (pelem[1]).wfract = 1.0-wfract1;
    (*pbulk).nr_elements=2 ;
    (*pbulk).massthickness = 100.0 ; /* lekker dik ! */
    (*pbulk).elementen=pelem;

   pmlayer->nr_layers=1;
   (*pmlayer).lagen[0]=pbulk;

    return(pmlayer);
}

multilayer *mktertbulk(char *elsym1, char *elsym2, char *elsym3,
                  double wfract1, double wfract2)
{
    layer *pbulk;
    elem *pelem;
    multilayer *pmlayer ;
   /*  create storage for layerdata and element dat */
 
   pbulk= (layer *)malloc(sizeof(layer));
   pelem=(elem *)malloc(3*sizeof(elem));
   pmlayer=(multilayer *)malloc(sizeof(multilayer));
   pmlayer->lagen=(layer **)malloc(sizeof(layer));
   strcpy((pelem[0]).name,elsym1) ;
   strcpy((pelem[1]).name,elsym2) ;
   strcpy((pelem[2]).name,elsym3) ;
   (pelem[0]).wfract = wfract1;
   (pelem[1]).wfract = wfract2;
   (pelem[2]).wfract = 1.0 -wfract1 -wfract2;
    (*pbulk).nr_elements=3 ;
    (*pbulk).massthickness = 100.0 ; /* lekker dik ! */
    (*pbulk).elementen=pelem;
 
   pmlayer->nr_layers=1;
   (*pmlayer).lagen[0]=pbulk;
    return(pmlayer);
}
