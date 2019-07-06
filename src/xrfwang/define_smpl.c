/******************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
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

void def_sample (multilayer *sample)
{
   FILE *fp;
   char filenaam[80];
   multicomplayer monster;
   printf("Give name of file with sample data ");
   gets(filenaam);
   fp=own_fopen(filenaam);
   bldcompos(fp,&monster);
   convcomp(&monster,sample);
}


