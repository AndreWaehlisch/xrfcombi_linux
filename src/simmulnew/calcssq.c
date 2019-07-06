/*
*   file calcssq.c
*
*   calculates sum of squared errors
*
*   from entries in 2 separate files
*   whose names should be given on the
*   command line
*
*   output to stdout
*
*   M.Bos
*   March 8 th 2001
*
*/
#include <stdio.h>
int main(int argc, char *argv[])
{
   int i=0;
   double result=0.0, number1, number2 ;

   FILE *fpin1, *fpin2 ;

   if ((fpin1=fopen(argv[1],"r"))==NULL)
      {
        printf("Cannot open %s \n", argv[1]);
        exit(1);
      }
   if (( fpin2=fopen(argv[2],"r"))==NULL)
      {
        printf("Cannot open %s \n", argv[2]);
        exit(1);
      }

  while ( fscanf(fpin1, "%lf", &number1) != EOF)
  {
      if (fscanf(fpin2, "%lf", &number2) == EOF)
      {
        printf("Second file has not enough numbers \n");
        exit(1);
      }
      /* debug 
      printf("%lf %lf \n",number1, number2); /* */
      result += (number1-number2)*(number1-number2);
      i++;
  } 
  printf("%lf", result);
  exit(0);
}

