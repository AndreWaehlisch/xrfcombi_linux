/*****************************************************
*
*	compdata.c
*	datstructure and access routines
*	for multilayer thin film samples
*	consisting of chemically defined 
* 	compounds
*
*	M.Bos
*	october 1997
*
*******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include "xrfluor.h"

extern FILE *own_fopen(char *);

/********************************************************
*  
* routine to build structure for sample data
* in terms of compounds 
*
*  inputs:
*  (1)  file pointer to file with
*  ascii representation of structure:
*  number of layers
*  per layer:
*  number of compounds, massthickness, compdata
*  compdata: compound symbol, weight fraction
*  (2) pointer to output variable of type multicomplayer
*
*  datstructures are mallocated from memory
*  and elements in the multicomplayer variable
*  are filled with adresses or data
*
********************************************************/

void bldcompos(FILE *fp, multicomplayer *alles)
{
	int nr_layers ;
   	int nr_comps ;
        int i,j, k, nr_atoms, cnt_at ;
        char *symbol_1;
        elem_cnt *el_pairs ;
        char regel[LINELENGTH];
        char *position;
        FILE *fpcomplist;
  	char sym_compound[MAX_CHAR_COMP] ;
	double weight_fraction, thickness;
        complayer **setlayers ;
        compound *reeks ;
    
        fscanf(fp, "%d", &nr_layers);
        /* debugging
        printf("nr of layers with compounds is %d\n",nr_layers); /* */
        setlayers = (complayer * * )malloc(nr_layers*sizeof(complayer));
        for (i=0; i < nr_layers; i++)
        {
  	  fscanf(fp,"%d",&nr_comps);
          fscanf(fp,"%lf",&thickness);
          /* debugging 
          printf("nr %d thickn %g fix %d \n",nr_comps, thickness,fix_thickn); /* */
	  reeks=(compound *)malloc(nr_comps*sizeof(compound));
          for (j=0; j < nr_comps; j++)
          {
             fpcomplist=own_fopen("./data/complist.dat");
             fscanf(fp,"%s %lf ",sym_compound, &weight_fraction);
             /* debugging 
             printf("gelezen is %s %g \n",sym_compound, weight_fraction); /* */
             reeks[j].wfract=weight_fraction;
             strcpy(reeks[j].name,sym_compound);

             /* now read from file complist to match compounds name */
 
             while ( (position=fgets(regel,LINELENGTH,fpcomplist)) != NULL)
             { 
               symbol_1=strtok(regel," ");   /* to symbol */
               if ( strcmp(symbol_1,sym_compound) ==0)
                break;
             }

             /* we now have the matching line or are past the end */

             if (position==NULL)
             {
               printf("unknown compound %s \n",sym_compound);
               exit(1);
             }

             position=strtok(NULL, " "); /* to density , throw away*/
             position=strtok(NULL," ");  /* AT or WT ? */
             if (strcmp(position, "AT") != 0)
             {
               printf("In complist.dat only AT entries for now\n"); 
               exit(1);
             }
             position=strtok(NULL," ");  /* to number of atoms */
             nr_atoms=atoi(position);
             el_pairs=(elem_cnt *) malloc(nr_atoms*sizeof(elem_cnt));
             reeks[j].nr_elems=nr_atoms;
             reeks[j].elementen=el_pairs;

             for (k=0; k < nr_atoms; k++)
             {
               position=strtok(NULL," ");   /* symbol of element */
               /* convert to lower case to match rest of system */
/*               if ((*(position+1) < 0x64) && ((*(position+1) != 0x20)) 
               *(position+1)= *(position+1)+0x20; */
               strcpy(el_pairs[k].name,position); /* name of element in place */
               position=strtok(NULL," \n");  /* atom count */
               cnt_at=atoi(position);
               el_pairs[k].count=cnt_at;
             }
           }
           setlayers[i]= (complayer *) malloc(sizeof(complayer));
           setlayers[i]->nr_compounds = nr_comps;
           setlayers[i]->comps=reeks;
           setlayers[i]->massthickness= thickness;
        }
        (*alles).nr_layers=nr_layers;
        (*alles).lagen=setlayers;
}
/*******
*  everything is now filled 
* here follow the access routines
*/

int get_nr_clay(multicomplayer *mclayer)
{
	return( (*mclayer).nr_layers);
}

int get_nr_comps(complayer *clayer)
{
    return(*clayer).nr_compounds;
}

complayer *get_clayer(multicomplayer *mclayer, int layer_nr)
{
      complayer *layer;
      layer = (*mclayer).lagen[layer_nr];
      return(layer);
}

char *get_comp_name(complayer *layer, int comp_nr)
{
     return ( (*layer).comps[comp_nr].name);
}

double get_cweightf(complayer *layer, int comp_nr)
{
    return( (*layer).comps[comp_nr].wfract);
}

int get_fixwght(complayer *layer, int comp_nr)
{
    return( (*layer).comps[comp_nr].wfixed);
}

double get_cmth(complayer *layer)
{
	return( (*layer).massthickness);
}

int get_fixmth(complayer *layer)
{
	return( (*layer).mthfixed);
}

/***************************
*
*  checkcomps()
*
*  for debugging purposes
* checks is access routines work
*  and that structures contain the required data
*
*************************************************/

void checkcomps(multicomplayer * mclayer)
{
	complayer *player ;
        int i, j ;
	int aant_comps ;
        int cnt_lay;
        cnt_lay=get_nr_clay(mclayer);
        printf("Number of layers with compounds is %d \n",cnt_lay);
        for (i=0; i < cnt_lay; i++) 
        {
	  player=get_clayer(mclayer, i);
          aant_comps=get_nr_comps(player);
          printf("Thickness layer %d (fix %d) is %g with %d compounds\n",
                  i, get_fixmth(player), get_cmth(player),aant_comps);
          for (j=0; j< aant_comps; j++)
          {
            printf("Compounds %s - fraction %g (fix %d) \n",
                    get_comp_name(player,j), get_cweightf(player,j),
                    get_fixwght(player,j));
           }
         }
}

 
