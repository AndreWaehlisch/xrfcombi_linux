/*
*
*       COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*	layerdata.c
*	datastructure and access routines
*	for multi-layer thin film samples
*
*	M.Bos
*	august 1997
*/
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "xrfluor.h"
  

/* ***************************************
*  routine to build a multilayer structure
*  
*  inputs :
*   (1) file pointer to file with 
*  ascii representation of structure:
*  number of layers
*  per layer:
*  number of elements, massthickness , el.data
*  el.data : Element symbol, weight fraction
*   (2) pointer to variable of type multilayer
*
*
*  datastructures are mallocated from memory
*  and elements in the multilayer variable are
*  filled with pertinent adresses or data
**************************************************/

void buildlay(FILE *fptest, multilayer *geheel)
{

  int aant_lagen;
  int  aant_elem ;
  int i, j;
  char sym_elem[4];
  double weight_fraction, dikte;
  layer **setlagen ;
  elem *reeks;
  layer *plaag;


   fscanf(fptest,"%d", &aant_lagen);
   /* debugging
   printf("aantal lagen is %d\n", aant_lagen); /* */
   setlagen=  ( layer * * ) malloc(aant_lagen*sizeof(layer));
   for (i=0; i< aant_lagen; i++)
   {
      fscanf(fptest,"%d",&aant_elem);
      fscanf(fptest,"%lf",&dikte);
      reeks=(elem *)malloc(aant_elem*sizeof(elem));
      for (j=0; j < aant_elem; j++)
      {
         fscanf(fptest,"%s %lf ",sym_elem,&weight_fraction);
         /* debugging
         printf("gelezen %s %g \n",sym_elem, weight_fraction); /* */
         reeks[j].wfract=weight_fraction; 
         strcpy(reeks[j].name,sym_elem);
      }
/* debugging
 
      for (j=0; j < aant_elem; j++)
      {
        
         printf("terug gelezen %s %g \n",reeks[j].name, reeks[j].wfract);
      }  /* */
      setlagen[i]= (layer *) malloc( sizeof(layer));
      setlagen[i]->nr_elements= aant_elem;
      setlagen[i]->elementen=reeks;
      setlagen[i]->massthickness=dikte;
   } 
   (*geheel).nr_layers=aant_lagen; 
   (*geheel).lagen=setlagen;
/* nu is de boel gevuld 
*  hoofdstructuur is "geheel" (type multilayer)
*  nu kijken of het er ook staat
*/

/* eerste item in die structure multilayer is
*  het aantal lagen
*/
   /* debugging
   printf("aantal lagen is %d \n", (*geheel).nr_layers); /* */
   for (i=0; i < (*geheel).nr_layers; i++)
   {
      plaag=(*geheel).lagen[i];
      /* debugging
      printf("pointer naar %d -e laag is %p\n",i,plaag); /* */
      aant_elem=plaag->nr_elements;
      /* debugging
      printf("aantal elementen in %d -e laag is %d\n",i,aant_elem); /*  */
      dikte=plaag->massthickness;
      /* debugging
      printf("dikte laag %d is %g \n",i, dikte);  
      for (j=0; j< aant_elem; j++)
      {
        printf("laag %d %s is %g  \n", i,(plaag->elementen[j]).name,
                (plaag->elementen[j]).wfract);
      } /* */
  
   }

}
/*********************************************************
*
*   routine used in debugging
*   prints all data in a multilayerstructure
*
*   inputs: pointer to a variable of type multilayer
*
***********************************************************/

void checkdata(multilayer *geheel)
{
  int i,j ;
  double dikte;
  int aant_lagen, aant_el ;
  layer *plaag ;
  elem *pelem;

  aant_lagen= (*geheel).nr_layers ;
  printf("aantal lagen is %d\n",aant_lagen);
  for (i=0; i < aant_lagen; i++)
  {
     plaag=(*geheel).lagen[i];
     aant_el=plaag->nr_elements;
     printf("laag %d heeft %d elementen\n", i, aant_el);
     dikte= plaag->massthickness;
     printf("massthickness laag %d is %g \n",i, dikte);
/*     printf("pelem is %p \n",pelem); */
     for (j=0; j < aant_el; j++)
     {
       pelem=plaag->elementen+j;
       printf("El is %s, fractie is %g \n",pelem->name, pelem->wfract);
     }
  }
}

/********************************************************
*
*   get_nr_layers(multilayer *) :
*   
*
*   returns the number of layers in an existing multilayer
*   structure
*
**********************************************************/

int get_nr_layers(multilayer *mlayer)
{
	return( (*mlayer).nr_layers);
}

/*********************************************************
*
*  get_nr_elements(layer *laag) :
*
*  returns the number of elements in a given layer
*
**********************************************************/

int get_nr_elements(layer *laag)
{
    return( (*laag).nr_elements);
}

/****************************************************
*
*  *get_layer()
*
*  returns the address of the layer data of a given
*  layer number within a multilayer
*
*  counting of the layers starts with zero at the top
*
********************************************************/
layer *get_layer(multilayer *mlaag, int laag_nr)
{
    layer *laag;
    laag= (*mlaag).lagen[laag_nr];
    return(laag);
}

/****************************************************
*
*   get_elem()
*
*   returns the address of the string of 
*   the chemical symbol
*   of a given element within a given layer
*
*****************************************************/
char *get_elem(layer *laag, int el_nr)
{
   return( (*laag).elementen[el_nr].name);
}

/****************************************************
*
*   get_nr_elem()
*
*   returns the number of the element with
*   given symbol  in a specific given layer
*   not in layer then return 200
*
********************************************************/

int get_nr_elem(layer *laag, char *symb)
{
    int i, tot_aant ;
    tot_aant = get_nr_elements(laag);
    i=0;
    while ( (i <tot_aant)&& (strcmp(get_elem(laag,i),symb))) i++ ;
    if (i==tot_aant)
    {
       printf("Element %s not in layer \n",symb);
       return(200);   /* element not present in layer */
    }
    return (i);
}


    
/***************************************************
*
*  get_weightfract()
*
*  returns the weight fraction of a given element
*  in a given layer
*
*******************************************************/
double get_weightfract(layer *laag, int el_nr)
{
   if (el_nr==200) return (0.0); else
   return( (*laag).elementen[el_nr].wfract);
}

/*********************************************************
*
* get_mthickness()
*
* returns the massthickness of a given layer
*
*********************************************************/
double get_mthickness(layer *laag)
{
   return( (*laag).massthickness);
}


/*******************************************************
*
*   check()
*
*   for debugging purposes to debug
*   the other acces routines
*
********************************************************/

void check(multilayer *mlaag)
{
    layer *plaag;
    int i,j ;
    int aantal_els;

    printf("Number of layers is %d\n", get_nr_layers(mlaag));
    for (i=0; i < get_nr_layers(mlaag); i++)
    {
        plaag=get_layer(mlaag,i);
        aantal_els=get_nr_elements(plaag);
        printf("Dikte laag %d is %g en hij telt %d elementen\n",
                i, get_mthickness(plaag),aantal_els);
        for (j=0; j < aantal_els; j++)
        {
          printf(" El %s  -- gehalte %g \n",get_elem(plaag,j),
                   get_weightfract(plaag,j));
        }
     }
}

