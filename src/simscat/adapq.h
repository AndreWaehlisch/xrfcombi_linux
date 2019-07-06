/***************************************************************************\
* Module : ADAPQ                                                            *
* Call   : Jan Bos                                                          *
* Date   : 09-04-1999                                     *
\***************************************************************************/
#ifndef _ADAPQ
#define _ADAPQ
#ifdef __cplusplus
extern "C" {                    /* fill in what we've left out */
#endif

typedef double (*extfun)(double,void *);
typedef double (*ext2fun)(double,double,void *);


double aquad(double xmin,double xmax,double Eps,extfun f,void *AddD);
double aquad2d(double xmin,double xmax,double ymin,double ymax,double Eps,ext2fun f,void *AddD);

#ifdef __cplusplus
 }
#endif
#endif
