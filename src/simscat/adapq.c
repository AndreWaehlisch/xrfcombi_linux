/***************************************************************************\
* Module      : ADAPQ                                                       *
* Call        : Jan                                                         *
* Purpose     : Adaptive integral method                                    *
* Date        : 07-02-94                                                    *
* Depends on  :                                                             *
* Remarks     :                                                             *
*                                                                           *
\***************************************************************************/
/*
#include "low.h"
#include "lowaddi.h"
#include "numadd.h"
#include "complx.h"
*/
#include <stdio.h>
#include <math.h>
#include "adapq.h"
#include "xrfluor.h"

/* #define M_PI (4.0*atan(1.0)) */

typedef struct {
 double (*ya)(double),(*yb)(double ),(*f)(double,double);
 } GArguments;

typedef struct {
 double ya, yb,(*f)(double,double);
 } GArguments2;

typedef struct {
 double x, (*f)(double , double);
 } FArguments;

typedef struct _F2D_1D {
 void *AddD2D;
 ext2fun f;
 double xb, xe, atx;
 double yb, ye, aty;
} F2D_1D;

static int    quadrule(int n, int n2,int m, int *maxrule, double *t, double *w,
                double *w1);
static double cheb(double x, double *a, int n);
static double fac(int n);
static double G (double x, GArguments  *gargs);
static double G2(double x, GArguments2 *gargs);
static double F (double y, FArguments  *fargs);
static double fac(int n);

static double fac(int n) {
 if (n==1) return(1.0);
 else return(1.0*n*fac(n-1));
 }

static double cheb(double x, double *a, int n)
       {
        int    r;
        double b0, b1, b2, twox=x+x;

        b1=0.0; b0=0.5*a[n];
        for (r=n-1; r>=1 ;r -= 1) {
            b2=b1; b1=b0;
            b0=twox*b1 - b2 +a[r];
        }
        return(x*b0-b1+0.5*a[0]);
       }

static int quadrule(int n, int n2, int m, int *maxrule, double *t, double *w, double *w1)
       {
 int    i;

 for (i=1; i<=n2 -1; i+=2) t[i*m]=cos(M_PI*i/n);
 for (i= *maxrule/2 + 2; i<=n ;i+=2)
     w1[i-1]=0.0,w1[i]=1.0/(i*i-1.0);
 for (i=1 ;i<=n2 ;i+=1)
     w[n2+i]=-4.0*cheb(t[i*m],w1,n)/n;
 *maxrule=n;
 return 0;
 }

/***************************************************************************\
* fun: int adapquad(extfun f,void *AddD,double a,double b, double *eps, double acc,double eta, int divmax, double *ans, double *err_r)*
*                                                                           *
* par:                                                                      *
*                                                                           *
* des: A doubly-adaptive Clenshaw-Curtis quadrature method J.Oliver.        *
*      The computer Journal Vol.15 No 2 pp.141-147 (No 42 BBV ref)          *
*                                                                           *
* com:                                                                      *
* ref:                                                                      *
* key: math                                                                 *
\***************************************************************************/
int adapquad(extfun f,void *AddD,double a,double b, double *eps, double acc,
             double eta, int divmax, double *ans, double *err_r) {
 int     i,j,m,mmax,n,n2,nmax, maxrule,order,dv,caution;
 double  c,cprev=0,e=0,eprev=0,fmax,fmin, h,hmin,I,intprev=0,k,k1,re,x,xa,
         xb,xc;
 double  ec[7], fs[65], xs[65], fx[1029], w1[1029],sigma[7][5], t[65],
         w[129];


 *ans=*err_r=0.0;
 if (a==b) goto COMPLETE;
 if (a>b) {
     c=b; b=a; a=c;
 }
 hmin=(b-a)/pow(2.0, divmax*1.0);
 if (acc<eta) acc=eta; acc *=16.0;
 x=4.0; xa= 64;
 for (i=1; i<=6; i++) {
     ec[i]=xa/((x*x -1.0)*(x*x -9.0));
     x+=x; xa+=xa;
 }

 /* critical values of k */

 sigma[1][0]=0.455; sigma[1][1]=0.272; sigma[1][2]=0.606;
 sigma[1][3]=0.811; sigma[1][4]=0.908;
 sigma[2][0]=0.55 ; sigma[2][1]=0.144; sigma[2][2]=0.257;
 sigma[2][3]=0.376; sigma[2][4]=0.511;
 sigma[3][0]=0.667; sigma[3][1]=0.243; sigma[3][2]=0.366;
 sigma[3][3]=0.449; sigma[3][4]=0.522;
 sigma[4][0]=0.78 ; sigma[4][1]=0.283; sigma[4][2]=0.468;
 sigma[4][3]=0.565; sigma[4][4]=0.624;
 sigma[5][0]=0.855; sigma[5][1]=0.290; sigma[5][2]=0.494;
 sigma[5][3]=0.634; sigma[5][4]=0.714;
 sigma[6][0]=-1.0 ; sigma[6][1]=0.292; sigma[6][2]=0.499;
 sigma[6][3]=0.644; sigma[6][4]=0.745;

 n=4; n2=2; nmax=128;
 m=mmax=nmax/n;
 t[0]=1.0; t[nmax/2]=0.0;
 maxrule=dv=0; w1[0]=-1.0;
 quadrule(n, n2, m, &maxrule, t, w, w1);
 xa=xc=a; xb=b;
 fx[0]=(*f)(b,AddD); fx[n]=(*f)(a,AddD);

NEXT:   /*printf("NEXT:\n");*/
 n=4; n2=2; m=mmax; order=1;
 caution=(xa<xc);
 h=xb-xa; k1=h/(b-xa);
 if (k1<0.1) k1=0.1;
 h*=0.5; j=1;
 fx[n2]=(*f)(xa+h,AddD);
 fmin=fmax=fx[0];
 if (fmax<fx[n]) fmax=fx[n];
 else {
   if (fmin>fx[n]) fmin=fx[n];
 }
 if (fmax<fx[n2]) fmax=fx[n2];
 else {
   if (fmin>fx[n2]) fmin=fx[n2];
 }

AGAIN:  /* printf("AGAIN:\n");/* */
 for (i=1;i<=n2-1; i+=j) {
     fx[i]=(*f)(xa+(1.0+t[i*m])*h,AddD);
     if (fmax<fx[i]) fmax=fx[i];
     else if (fmin>fx[i]) fmin=fx[i];
     fx[n-i]=(*f)(xa+(1.0-t[i*m])*h,AddD);
     if (fmax<fx[n-i]) fmax=fx[n-i];
     else if (fmin>fx[n-i]) fmin=fx[n-i];
 }
/* if (fabs(fmax)>fabs(fmin)) re=acc*fabs(fmax);
 else re=acc*fabs(fmin);*/
 re=(fabs(fmax)>fabs(fmin))?acc*fabs(fmax):acc*fabs(fmin);
/* if (n==4) j= 4 ;
 else j=6;*/
 j=(n==4)?4:6;
 k=0.0;
 c=fabs(cheb(-1.0, fx,n))/n;
 for (i=2; i<=j; i+=2) {
/*     if (i<=n2) x= -t[i*m];
     else x=t[(n-i)*m];*/
     x=(i<=n2)?-t[i*m]:t[(n-i)*m];
     cprev=c; c=fabs(cheb(x, fx,n))/n2;
     if (c>re) {
      if (k<cprev/c) k=cprev/c;
      }
     else {
      if (cprev > re) k=1.0;
      }
 } /* nieuwe einde loop */
 if (k>sigma[order][4]) {
    if (n==4) goto SPLIT;
    else goto EVAL;
 }
 if (n==4) cprev=c;
 else {
   if (cprev<re) cprev=k*c;
 }
/* } was eind loop */
 e=h*cprev*ec[order]*k*k*k;
 for (i=0; (k>sigma[order][i]); i+=1) e*=2.0;
 re*=h;

EVAL:       /* printf("EVAL\n"); /* */
 I=w1[n]*(fx[0]+fx[n])+w[n]*fx[n2];
 for (i=1; i<=n2-1; i+=1) I+=w[n2+i]*(fx[i]+fx[n-i]);
 I*=h;
 if (n==4) goto TEST;
 c=fabs(I-intprev);
 if (c>eprev) {
     caution=1;
     if (xc<xb) xc=xb;
 }
 else caution=0;
 if ((k>sigma[order][4])||caution) e=c;
 if (e>c) e=c;

TEST:      /* printf("TEST:\n");/* */
 if ((e<re)||(e<=k1*(*eps))) goto UPDATE;
 if (k>sigma[order][0]) goto SPLIT; else goto DOUBLE_tag;

UPDATE:    /* printf("UPDATE\n"); /* */
 if ((n==4)&&(caution||((xa==a)&&(dv==0)))) goto DOUBLE_tag;
 if (e<re) e=re;
 *err_r+=e; *eps-=e;
 if (*eps<0.1* (*err_r)) *eps=0.1* (*err_r);
 *ans+=I;
 if (dv==0) goto COMPLETE;
 dv--;
 xa=xb; xb=xs[dv];
 fx[4]=fx[0]; fx[0]=fs[dv];
 goto NEXT;

DOUBLE_tag:  /*   printf("DOUBLE:\n");/* */
 for (i=n ; i>=1 ; i -= 1) fx[2*i]=fx[i];
 n2=n; n+=n; m/=2;
 order++;
 eprev=e; intprev=I;
 if (eprev<re) eprev=re;
 if (n>maxrule) quadrule(n, n2, m, &maxrule, t, w, w1);
 j=2;
 goto AGAIN;

SPLIT:    /*  printf("SPLIT:\n"); /* */
 e=2.0*h*(fmax-fmin);
 if((e<re)||(e<=k1*(*eps))) {
     I=h*(fmax+fmin);
     goto UPDATE;
 }
 if (h<hmin) return(1);
 xs[dv]=xb; xb=xa+h;
 fs[dv]=fx[0]; fx[0]=fx[n2]; fx[4]=fx[n];
 dv++;
 goto NEXT;

COMPLETE:   /* printf("COMPLETE:\n");/* */
 return(0);
 }

/***************************************************************************\
* fun: double aquad(double xmin,double xmax,double Eps,extfun f,void *AddD) *
*                                                                           *
* par:                                                                      *
*                                                                           *
* des: Adaptive integral method, basically an adaption on quad.             *
*                                                                           *
*                                                                           *
* com:                                                                      *
* ref:                                                                      *
* key: math                                                                 *
\***************************************************************************/
static double __aquad(double xmin,double xmax,double Eps,
                      extfun f,void *AddD,int DivSize) {
 double d,eta=eta=pow(2.0,-64.0),err_r=0.0;
 int re;
 double f1,f2;
 f1=f(xmin,AddD);
 f2=f(xmax,AddD);
 if((fabs(f1-f2)<1e-15) && DivSize) {
  double xmi=(xmin+xmax)/2;
  return __aquad(xmin,xmi,Eps,f,AddD,DivSize/2)+
         __aquad(xmi,xmax,Eps,f,AddD,DivSize/2);
  }
 re=adapquad(f,AddD,xmin,xmax,&Eps,0.0,eta,512,&d,&err_r);
 if (re) printf("Error in adapquad!\n");
 return re?0.0:d;
 }

double aquad(double xmin,double xmax,double Eps,extfun f,void *AddD) {
 return __aquad(xmin,xmax,Eps,f,AddD,128);
 }

static double F(double y, FArguments *fargs) {
 return(fargs->f(fargs->x,y));
 }

static double G(double x, GArguments *gargs) {
 double yai, ybi,ans, err_r, EPS=1.0e-14;
 FArguments fargs;
 yai=gargs->ya(x);
 ybi=gargs->yb(x);
 fargs.x=x;
 fargs.f=gargs->f;
 adapquad((extfun)F,&fargs,yai,ybi,&EPS,10e-14,pow(2.0,-64.0),64,&ans,&err_r);
 return(ans);
 }

static double G2(double x, GArguments2 *gargs) {
 double yai, ybi,ans, err_r, EPS=1.0e-14;
 FArguments fargs;
 yai=gargs->ya;
 ybi=gargs->yb;
 fargs.x=x;
 fargs.f=gargs->f;
 adapquad((extfun)F,&fargs,yai,ybi,&EPS,10e-14,pow(2.0,-64.0), 64, &ans, &err_r);
 return(ans);
 }

double var_int2d(double xa, double xb, double ( *lb)(double),
 double ( *ub)(double), double (*f)(double,double),double EPS) {
 double ans, err_r;
 GArguments gargs;
 gargs.ya=lb;
 gargs.yb=ub;
 gargs.f =f;
 adapquad((extfun)G, &gargs, xa, xb, &EPS, 10e-14,
               pow(2.0,-64.0), 64, &ans, &err_r);
 return ans;
 }

double int2d(double xa, double xb, double ya, double  yb, double (*f)(double,double), double EPS) {
 double ans, err_r;
 GArguments2 gargs;
 gargs.ya=ya;
 gargs.yb=yb;
 gargs.f =f;
 adapquad((extfun)G2, &gargs, xa, xb, &EPS, 10e-14,
               pow(2.0,-64.0), 64, &ans, &err_r);
 return ans;
 }

static double F_x0_y(double y, void *Add2D1D) {
 F2D_1D *AddD= (F2D_1D *) Add2D1D;
 return(AddD->f(AddD->atx,y,AddD->AddD2D));
}

static double F_x_y0(double x, void *Add2D1D) {
 F2D_1D *AddD= (F2D_1D *) Add2D1D;
 return(AddD->f(x,AddD->aty,AddD->AddD2D));
}

static double Inty_fdy(double x, void *Add2D1D) {
 double  Eps=1.0e-12;
 double d,eta=eta=pow(2.0,-64.0),err_r=0.0;
 int re;
 F2D_1D *AddD= (F2D_1D *) Add2D1D;
 AddD->atx=x;
 re=adapquad((extfun) F_x0_y ,(void *)AddD,AddD->yb,AddD->ye,&Eps,1e-14,eta,512,&d,&err_r);
/* ans=aquad(AddD->yb,AddD->ye, EPS,(extfun) F_x0_y,(void *) AddD) ;
 return(ans);*/
 if (re) printf("Error in adapquad!\n");
 return re?0.0:d;

 }

double aquad2d(double xa, double xb, double ya, double  yb, double Eps, ext2fun f, void *AddD) {
 F2D_1D AddD1D;
 double d,eta=eta=pow(2.0,-64.0),err_r=0.0;
 int re;

 AddD1D.AddD2D=AddD;
 AddD1D.f=f;
 AddD1D.xb=xa;
 AddD1D.xe=xb;
 AddD1D.yb=ya;
 AddD1D.ye=yb;

/* ans=aquad(xa,xb,Eps,(extfun) Inty_fdy,(void *)&AddD1D) ;
 return(ans);*/
 re=adapquad((extfun) Inty_fdy,(void *)&AddD1D,xa,xb,&Eps,1e-14,eta,512,&d,&err_r);
 if (re) printf("Error in adapquad!\n");
 return re?0.0:d;

}


