/*
*  file fibonacci.c
*
*
*  calculates  offset in interval
*  using Fibonacci series
*   1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233
*
*  input values low, high, tolerance as commandline args
*
*  M.Bos
*  March 2001
*
*/
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  int i, fibo[20], nr_exps, n  ;
  double x1, x2, l1, l2 ;
  double low, high, tolerance, breuk;

  if ( argc != 4)
  {
    printf("usage: %s low high tolerance stepnr\n");
    exit(1);
  }
  
  low= atof(argv[1]);
  high=atof(argv[2]);
  tolerance=atof(argv[3]);
  /* printf("low %lf high %lf tol %lf stepnr %d \n", low, high,
         tolerance); */
  fibo[0]=1;
  fibo[1]=1;
  i=1;

  /* preliminary calcs for start position */

    breuk= (high-low)/tolerance;

    while ( breuk > fibo[i++])
    {
       fibo[i]=fibo[i-1]+fibo[i-2] ;
    /*   printf("fibo[%d] is %d \n", i , fibo[i]); */
    }
  n=i-1 ;
  l1= ((double)fibo[n-2]/(double)fibo[n])*(high-low);
  printf("%lf\n",l1); 
  exit(0);
}

