/*
this is the main program to transform coordinates. it sends points
from the transformed coordinates (in the normal form or central
manifold) to the initial ones and viceversa. it is assumed that the
changes have been obtained (and written into several files after
realification) by the others programs of the package. the files with
the expansions of the changes of variables have to be in the directory
contained in the variable DATA (see makefile). the program needs three
parameters in the command line: the name of the file with the initial
points (input, assumed to be in the directory pointed by DATA), the
file where the transformed coordinates will be written (output, to be
created in the directory pointed by DATA), and an integer value equal
to 1 or -1. this value determines the direction of the transformation
(if 1, the initial points are assumed to be in normal form coordinates
and the final one in initial coordinates; and viceversa for -1). if
the variable DATA is not defined in the compiler command line (see
makefile), all the files are assumed to be in the current directory.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arit-r.h"

/*
now we define the complex type. although it is not necessary in this
program, we need the definition of complex to avoid compiler complains
when processing files io6s.h and io6p.h, since the type complex is
used there in the prototypes of some functions.
*/

typedef struct {
renum re;
renum im;
} complex;

#include "mp4s.h"
#include "mp6s.h"
#include "mp6p.h"
#include "io4s.h"
#include "io6s.h"
#include "io6p.h"
#include "iol.h"
#include "uneix.h"
#include "erpe6.h"
#include "repo.h"
#include "litr5.h"

/*
next lines are to ensure that the variable DATA is defined. it
must contain the directories where data files are going to be
stored/found. it is usually defined in the command line (see
makefile).
*/

#if !defined(DATA)
#define DATA ""
#endif

#include "files-nf.h"

int main(int narg, char *arg[])
{
   void litr123(renum a[6][6], integer l, renum mu, renum gam,
             renum x[6], renum y[6], integer dir);
   /*void litr123_(renum a[4][4], integer l, renum mu, renum gam,
             renum x[4], renum y[4], integer dir); */
   void explain(char *s);
   renum a[6][6],y[6],**c,**pi,**pf,e,mu,gam;
   //renum a[4][4],y[4],**c,**pi,**pf,e,mu,gam;
   double *erro;
   integer i,j,li,n,np,m,ms,mp,dir;
   char *cc[6],*ipo,*fpo,*cvl,*ctl; /* filenames involved */
   FILE *f;
   if (narg != 4) {explain(DATA); exit(1);}
   ipo=uneix(DATA,arg[1]);
   fpo=uneix(DATA,arg[2]);
   dir=atoi(arg[3]);
   if ((dir != 1) && (dir != -1)) {explain(DATA); exit(1);}
/*
   next lines are to join the path where the files have to be found/stored
   (given by DATA) and the name of the files (defined in files-nf.h)
*/
   if (dir == 1)
      {
         cc[0]=uneix(DATA,CD1);
         cc[1]=uneix(DATA,CD2);
         cc[2]=uneix(DATA,CD3);
         cc[3]=uneix(DATA,CD4);
         cc[4]=uneix(DATA,CD5);
         cc[5]=uneix(DATA,CD6);
      }
      else
      {
         cc[0]=uneix(DATA,CI1);
         cc[1]=uneix(DATA,CI2);
         cc[2]=uneix(DATA,CI3);
         cc[3]=uneix(DATA,CI4);
         cc[4]=uneix(DATA,CI5);
         cc[5]=uneix(DATA,CI6);
      }
   cvl=uneix(DATA,CVL);
   ctl=uneix(DATA,CTL);
   rcvl(a,cvl);
   rctl123(ctl,&li,&n,&mu,&gam);
   printf("equilibrium point: %d\n",li);
   printf("mass parameter: %24.16e\n",mu);
   //puts("degree of the change of variables?");
   //scanf("%d",&n); /* we overwrite the variable n read by rctl123 */
/*
   next lines are to pass the mass parameter to the program that will
   perform the numerical integration for the rtbp
*/
   if (dir == 1)
      {
         ctl=uneix(DATA,"rtbp.mu");
         f=fopen(ctl,"w");
         if (f == NULL) {printf("main-tcnf: cannot open %s\n",ctl); exit(1);}
         fprintf(f,"%25.16e\n",mu);
         fclose(f);
      }
   //imp4s(n);
   imp6s(n);
   imp6p(n);
   c=(renum**)malloc((n+1)*sizeof(renum*));
   if (c == NULL) {puts("main-tcnf: out of memory (1)"); exit(1);}
   for (i=0; i<=n; i++)
   {
      //m=ntph4s(i);
      ms=ntph6s(i);
      mp=ntph6p(i);
      m = ((ms > mp) ? ms : mp);
      c[i]=(renum*)malloc(m*sizeof(renum));
      if (c[i] == NULL) {puts("main: out of memory 2"); exit(1);}
   }
   pi=repo(ipo,6,&np); /* to read the initial points */
   printf("number of read points: %d\n",np);
   erro=(double*)malloc(np*sizeof(double));
   for(j=0;j<np;j++) erro[j]=0.;
   pf=(renum**)malloc(np*sizeof(renum*));
   if (pf == NULL) {puts("main-tcnf: out of memory (2)"); exit(1);}
   for (j=0; j<np; j++)
   {
      pf[j]=(renum*)malloc(6*sizeof(renum));
      if (pf[j] == NULL) {puts("main-tcnf: out of memory (3)"); exit(1);}
   }
   if (dir == -1)
      {
         for (j=0; j<np; j++)
         {
            litr123(a,li,mu,gam,pi[j],y,dir);
            for (i=0; i<6; i++) pi[j][i]=y[i];
         }
      }
   for (i=0; i<4; i++)
   {
      //rreb6s(c,1,n,cc[i]);
      rrea6s(c,1,n,cc[i]);
      for (j=0; j<np; j++)
      {
         ini_erpe6(pi[j],n);
         //printf("cheguei aqui?\n");
         pf[j][i]=erpe6s(c,n,&e); //printf("%03d %03d %25.16e\n",i,j,e); 
         erro[j]+=fabs(e);
         end_erpe6();
      }
   }
   for (i=4; i<6; i++)
   {
      //rreb6p(c,1,n,cc[i]);
      rrea6p(c,1,n,cc[i]);
      for (j=0; j<np; j++)
      {
         ini_erpe6(pi[j],n);
         pf[j][i]=erpe6p(c,n,&e);
         end_erpe6();
      }
   }
   if (dir == 1)
      {
         for (j=0; j<np; j++)
         {
            litr123(a,li,mu,gam,pf[j],y,dir);
            for (i=0; i<6; i++) pf[j][i]=y[i];
         }
      }
   //for(j=0;j<np;j++) for(i=0;i<6;i++) pf[j][i]=pi[j][i];
   f=fopen(fpo,"w");
   if (f == NULL) {printf("main-tcnf: cannot open %s\n",fpo); exit(1);}
   fprintf(f,"%d %d\n",np,6);
   for (j=0; j<np; j++)
   {
      for (i=0; i<6; i++) fprintf(f," %24.16e",pf[j][i]);
      fprintf(f," %24.16e",erro[j]);
      fprintf(f,"\n");
   }
   fclose(f);
/*
 that's all.
*/
   for (i=0; i<np; i++) free(pf[i]);
   free(pf);
   for (i=0; i<np; i++) free(pi[i]);
   free(pi);
   for (i=0; i<=n; i++) free(c[i]);
   free(c);
   //amp4s();
   amp6p();
   amp6s();
   for (i=0; i<6; i++) free(cc[i]); /* it was allocated inside uneix */
   return(0);
}
void explain(char *s)
{
   puts("\nyou must provide three values in the command line:");
   puts("the names of two files and and a flag.");
   puts("the first file contains the initial points and the second");
   puts("one is the file where the transformed points are going");
   puts("to be stored");
   puts("the flag must be 1 or -1. 1 means that we go from normal");
   puts("form coordinates to initial ones, and -1 means the opposite");
   puts("transformation.\n");
   if (strlen(s) == 0)
      {
         puts("both files are assumed to be in the actual directory");
      }
      else
      {
         printf("both files are assumed to be in %s\n\n",s);
      }
   return;
}

void litr123(renum a[6][6], integer l, renum mu, renum gam,
             renum x[6], renum y[6], integer dir)
/*
this is to transform points coming from the nonlinear change
corresponding to the central manifold computation, into synodical
coordinates and vice-versa.

parameters:
a:   change used to diagonalize h2 (input).
x:   initial point (input).
y:   point in synodical coordinates (output).
dir: direction of the transformation. 1 means that the change is given
     by matrix a, -1 means that it is given by the inverse of a. as a
     is a symplectic matrix, we will compute its inverse using this
     property.
*/
{
   int i;
   renum c,w[6],z[6],v[6];

   if(dir==1){
      for (i=0; i<6; i++)
         w[i]=a[i][0]*x[0]+a[i][1]*x[1]+a[i][2]*x[2]
             +a[i][3]*x[3]+a[i][4]*x[4]+a[i][5]*x[5];

      z[0]=w[0];   z[1]=w[2]+w[1];
      z[2]=w[2];   z[3]=-w[0]+w[3];
      z[4]=w[4];   z[5]=w[5];

      switch (l)
      {
         case 1:
            c=-1+gam;
            w[0]=-gam*z[0]+mu+c;    w[1]=-gam*z[1];
            w[2]=-gam*z[2];         w[3]=-gam*z[3];
            w[4]=gam*z[4];          w[5]=gam*z[5];
            break;
         case 2:
            c=-1-gam;
            w[0]=-gam*z[0]+mu+c;    w[1]=-gam*z[1];
            w[2]=-gam*z[2];         w[3]=-gam*z[3];
            w[4]=gam*z[4];          w[5]=gam*z[5];
            break;
         case 3:
            c=gam;
            w[0]=gam*z[0]+mu+c;     w[1]=gam*z[1];
            w[2]=gam*z[2];          w[3]=gam*z[3];
            w[4]=gam*z[4];          w[5]=gam*z[5];
            break;
         default:
            printf("tic error: li=%d\n",l);
            exit(1);
      }

      y[0]=w[0];
      y[1]=w[1]-w[2];
      y[2]=w[2];
      y[3]=w[3]+w[0];
      y[4]=w[4];
      y[5]=w[5];
   }else if(dir==-1){

      /*z[0]=x[0];   z[1]=x[1]-x[2];
      z[2]=x[2];   z[3]=x[0]+x[3];
      z[4]=x[4];   z[5]=x[5];

      switch (l)
      {
         case 1:
            c=-1+gam;
            w[0]=-(z[0]-mu-c)/gam;  w[1]=-z[1]/gam;
            w[2]=-z[2]/gam;         w[3]=-z[3]/gam;
            w[4]=z[4]/gam;          w[5]=z[5]/gam;
            break;
         case 2:
            c=-1-gam;
            w[0]=-(z[0]-mu-c)/gam;  w[1]=-z[1]/gam;
            w[2]=-z[2]/gam;         w[3]=-z[3]/gam;
            w[4]=z[4]/gam;          w[5]=z[5]/gam;
            break;
         case 3:
            c=gam;
            w[0]=(z[0]-mu-c)/gam;   w[1]=z[1]/gam;
            w[2]=z[2]/gam;          w[3]=z[3]/gam;
            w[4]=z[4]/gam;          w[5]=z[5]/gam;
            break;
         default:
            printf("tic error: li=%d\n",l);
            exit(1);
      }

      y[0]=w[0];
      y[1]=w[1]+w[2];
      y[2]=w[2];
      y[3]=w[3]-w[0];
      y[4]=w[4];
      y[5]=w[5];*/

      switch(l){
         case 1:
            c=-1+gam;
            y[0]=-(x[0]-mu-c)/gam;  y[1]=-x[1]/gam;
            y[2]=-x[2]/gam;         y[3]=-(x[3]-mu-c)/gam;
            y[4]=x[4]/gam;          y[5]=x[5]/gam;
            break;
         case 2:
            c=-1-gam;
            y[0]=-(x[0]-mu-c)/gam;  y[1]=-x[1]/gam;
            y[2]=-x[2]/gam;         y[3]=-(x[3]-mu-c)/gam;
            y[4]=x[4]/gam;          y[5]=x[5]/gam;
	    //printf("%25.16e%25.16e%25.16e\n%25.16e%25.16e%25.16e\n",y[0],y[1],y[2],y[3],y[4],y[5]); exit(1);
            break;
         case 3:
            c=gam;
            y[0]=(x[0]-mu-c)/gam;   y[1]=x[1]/gam;
            y[2]=x[2]/gam;          y[3]=(x[3]-mu-c)/gam;
            y[4]=x[4]/gam;          y[5]=x[5]/gam;
            break;
         default:
            printf("tic error: li=%d\n",l);
            exit(1);
      }

         v[0]=y[1]; v[1]=-y[0];
         v[2]=y[3]; v[3]=-y[2];
         v[4]=y[5]; v[5]=-y[4];

         for (i=0; i<6; i++)
            y[i]=a[0][i]*v[0]+a[1][i]*v[1]+a[2][i]*v[2]
                +a[3][i]*v[3]+a[4][i]*v[4]+a[5][i]*v[5];

         v[0]=y[1]; v[1]=-y[0];
         v[2]=y[3]; v[3]=-y[2];
         v[4]=y[5]; v[5]=-y[4];

         for (i=0; i<6; i++) y[i]=-v[i];
   }else{
      printf("tic error: dir=%d\n",dir);
      exit(1);
   }

   return;
}
void litr123_(renum a[4][4], integer l, renum mu, renum gam,
              renum x[4], renum y[4], integer dir)
/*
this is to transform points coming from the nonlinear change
corresponding to the central manifold computation, into synodical
coordinates and vice-versa.

parameters:
a:   change used to diagonalize h2 (input).
x:   initial point (input).
y:   point in synodical coordinates (output).
dir: direction of the transformation. 1 means that the change is given
     by matrix a, -1 means that it is given by the inverse of a. as a
     is a symplectic matrix, we will compute its inverse using this
     property.
*/
{
   int i;
   renum c,w[4],z[4],v[4];

   if(dir==1){
      for (i=0; i<4; i++)
         w[i]=a[i][0]*x[0]+a[i][1]*x[1]
             +a[i][2]*x[2]+a[i][3]*x[3];

      z[0]=w[0];   z[1]=w[2]+w[1];
      z[2]=w[2];   z[3]=-w[0]+w[3];

      switch (l)
      {
         case 1:
            c=-1+gam;
            w[0]=-gam*z[0]+mu+c;    w[1]=-gam*z[1];
            w[2]=-gam*z[2];         w[3]=-gam*z[3];
            break;
         case 2:
            c=-1-gam;
            w[0]=-gam*z[0]+mu+c;    w[1]=-gam*z[1];
            w[2]=-gam*z[2];         w[3]=-gam*z[3];
            break;
         case 3:
            c=gam;
            w[0]=gam*z[0]+mu+c;     w[1]=gam*z[1];
            w[2]=gam*z[2];          w[3]=gam*z[3];
            break;
         default:
            printf("tic error: li=%d\n",l);
            exit(1);
      }

      y[0]=w[0];
      y[1]=w[1]-w[2];
      y[2]=w[2];
      y[3]=w[3]+w[0];
   }else if(dir==-1){

      switch(l){
         case 1:
            c=-1+gam;
            y[0]=-(x[0]-mu-c)/gam;  y[1]=-x[1]/gam;
            y[2]=-x[2]/gam;         y[3]=-(x[3]-mu-c)/gam;
            break;
         case 2:
            c=-1-gam;
            y[0]=-(x[0]-mu-c)/gam;  y[1]=-x[1]/gam;
            y[2]=-x[2]/gam;         y[3]=-(x[3]-mu-c)/gam;
            break;
         case 3:
            c=gam;
            y[0]=(x[0]-mu-c)/gam;   y[1]=x[1]/gam;
            y[2]=x[2]/gam;          y[3]=(x[3]-mu-c)/gam;
            break;
         default:
            printf("tic error: li=%d\n",l);
            exit(1);
      }

         v[0]=y[1]; v[1]=-y[0];
         v[2]=y[3]; v[3]=-y[2];

         for (i=0; i<4; i++)
            y[i]=a[0][i]*v[0]+a[1][i]*v[1]
                +a[2][i]*v[2]+a[3][i]*v[3];

         v[0]=y[1]; v[1]=-y[0];
         v[2]=y[3]; v[3]=-y[2];

         for (i=0; i<4; i++) y[i]=-v[i];
   }else{
      printf("tic error: dir=%d\n",dir);
      exit(1);
   }

   return;
}
