/******************************************
*
*   COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*   define_smpl.c
*
*   reads datafile  from file
*   and stores sample data in memory allocated
*   structure sample is filled with pointers
*   to the memory allocated
*
*   M.Bos
*   sept 1997
*
*******************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

extern FILE *own_fopen(char *);

void def_sample(FILE *fp, multilayer *sample, multicomplayer *monster)
{
   bldcompos(fp,monster);
   /* debugging
   printf("na bldcompos in define_smpl \n"); /* */
   convcomp(monster,sample);
}


