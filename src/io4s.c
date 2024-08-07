/*
this file contains routines that read polynomials handled by mp4s from
a binary or ascii file. note that the definition of 'complex' is done
'by hand', instead of including the file arit-c.h. so, a complex is a
structure containing two doubles.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp4s.h"
#include "mp6s.h"
#include "msgs.h"

typedef struct {
   double re;
   double im;
} complex;

#define ERR -1 /* value returned by fprintf on failure */

void wea4s(complex **h, integer ni, integer nf, char *nom)
/*
this is to write, in ascii format, the expansion h in the file nom.
the file is also created by this routine. if the file already exists,
it is deleted first.

parameters:
h:   expansion to be written (input).
ni:  degree to start writting the hamiltonian.
nf:  last degree to be written.
nom: string with the name of the file to write h (input).
*/
{
   integer i,j,k[4],l,m,c;
   char *fmt;
   FILE *f;
   fmt="%2d %2d %2d %2d %25.16e %25.16e\n";
   printf("wea4s: creating ascii file %s...\n",nom);
   f=fopen(nom,"w");
   if (f == NULL) {printf("wea4s: can't open file %s\n",nom); exit(1);}
   l=0;
   for (i=ni; i<=nf; i++)
   {
      m=ntph4s(i);
      for (j=0; j<m; j++)
      {
         llex4s(j,k,i);
         c=fprintf(f,fmt,k[0],k[1],k[2],k[3],h[i][j].re,h[i][j].im);
         if (c == ERR)
            {
               fclose(f);
               iomsg_ea2("wea4s",l); /* this will stop the program */
            }
         ++l;
      }
   }
   printf("total number of (written) monomials: %d\n",l);
   printf("closing %s\n",nom);
   fclose(f);
   return;
}
void rea4s(complex **h, integer ni, integer nf, char *nom)
/*
this is to read a complex expansion stored in the ascii file nom.
warning: before reading the file, the routine will fill with zeros
the expansion c, from degree 0 (0, not ni) to nf. if you don't
want this, you should change the marked loop below.

parameters:
h:   read expansion (output).
n:   degree of the expansion h.
nom: string with the name of the file to read (input).
*/
{
   double u,v;
   integer i,j,k[4],m;
   char  *fmt;
   FILE *f;
   fmt="%d %d %d %d %le %le";
   printf("rea4s: opening ascii file %s...\n",nom);
   f=fopen(nom,"r");
   if (f == NULL) {printf("rea4s: can't open file %s\n",nom); exit(1);}
   for (i=0; i<=nf; i++) /* this is the loop that zeroes h */
   {
      m=ntph4s(i);
      for (j=0; j<m; j++) {h[i][j].re=0.e0; h[i][j].im=0.e0;}
   }
   m=0;
   while (fscanf(f,fmt,k,k+1,k+2,k+3,&u,&v) != EOF)
   {
      i=k[0]+k[1]+k[2]+k[3];
      if ((i >= ni) && (i <= nf))
         {
            j=exll4s(k,i);
            h[i][j].re=u;
            h[i][j].im=v;
            ++m;
         }
   }
   printf("rea4s: total number of (read) monomials: %d\n",m);
   printf("rea4s: closing %s\n",nom);
   fclose(f);
   return;
}
void web4s(complex **h, integer ni, integer nf, char *nom)
/*
this is to write, in binary format, the expansion h in the file nom.
the file is also created by this routine. if the file already exists,
it is deleted first.

parameters:
h:   expansion to be written (input).
ni:  degree to start writting the hamiltonian.
nf:  last degree to be written.
nom: string with the name of the file to write h (input).
*/
{
   integer i,j,l,m,inf[4];
   FILE *f;
   printf("web4s: creating binary file %s...\n",nom);
   f=fopen(nom,"wb");
   if (f == NULL) {printf("web4s: can't open file %s\n",nom); exit(1);}
   inf[0]=4;  /* number of variables */
   inf[1]=2;  /* kind of simmetry */
   inf[2]=ni; /* initial degree */
   inf[3]=nf; /* final degree */
   j=fwrite(inf,sizeof(integer),4,f);
   if (j < 4) iomsg_eb2("web4s",nom); /* this will stop the program */
   l=0;
   for (i=ni; i<=nf; i++)
   {
      m=ntph4s(i);
      j=fwrite(h[i],sizeof(complex),m,f);
      if (j < m) iomsg_eb2("web4s",nom); /* this will stop the program */
      l += j;
   }
   printf("web4s: total number of (written) monomials: %d\n",l);
   printf("web4s: closing %s\n",nom);
   fclose(f);
   return;
}
void reb4s(complex **h, integer ni, integer nf, char *nom)
/*
this is to read an expansion contained in the binary file nom.

parameters:
h:   read series (output).
ni:  first degree to be read.
nf:  last degree to be read.
nom: name of the binary file containing the expansion (input).
*/
{
   integer i,j,k,l,m,inf[4];
   complex w;
   FILE *f;
   printf("reb4s: opening binary file %s...\n",nom);
   f=fopen(nom,"rb");
   if (f == NULL) {printf("reb4s: can't open file %s\n",nom); exit(1);}
   j=fread(inf,sizeof(integer),4,f); /* this is the head of the file */
   if (j < 4) {puts("reb4s error 1."); exit(1);}
   if ((inf[0] != 4) || (inf[1] != 2)) iomsg_eb02("reb4s",nom,inf[0],inf[1]);
   iomsg_eb01("reb4s",nom,ni,nf,inf[2],inf[3]);
/*
   the procedure used to read the file is the following: first we set
   the array c (for degrees between ni and nf) to zero and then we read
   the monomials of the file whose degree belong to the interval [ni,nf].
   in this way, all the monomials that are not in the file are considered
   to be zero.
*/
   for (i=ni; i<=nf; i++)
   {
      m=ntph4s(i);
      for (j=0; j<m; j++) {h[i][j].re=0; h[i][j].im=0;}
   }
   l=0;
   for (i=inf[2]; i<=inf[3]; i++)
   {
      if (i > nf) break;
      if (i >= ni)
         {
            m=ntph4s(i);
            j=fread(h[i],sizeof(complex),m,f);
            if (j < m) {puts("reb4s: reading error 1"); exit(1);}
            l += j;
         }
         else
         {
            m=ntph4s(i);
            for (j=0; j<m; j++)
            {
               k=fread(&w,sizeof(complex),m,f);
               if (k < 1) {puts("reb4s: reading error 2"); exit(1);}
            }
         }
   }
   printf("reb4s: total number of (read) monomials: %d\n",l);
   printf("reb4s: closing %s\n",nom);
   fclose(f);
   return;
}
void wcm4s(complex **h, integer ni, integer nf, char *nom, double tol)
/*
this is to write, in ascii format, the realified expansion for the
central manifold in the file nom. only real parts bigger than tol are
written. if some imaginary part is bigger than tol, a warning message
is issued. the file is also created by this routine. if the file
already exists, the previous one is deleted.

parameters:
h:   expansion to be written (input).
ni:  degree to start writting the hamiltonian.
nf:  last degree to be written.
nom: string with the name of the file to write h (input).
tol: threshold to decide what is zero and what is not.
*/
{
   integer i,j,k[4],l,m,e;
   double a,b,c,d;
   char *fm1,*fm2;
   FILE *f;
   fm1="wcm4s warning: %2d %2d %2d %2d %25.16e %25.16e\n";
   fm2="%2d %2d %2d %2d %25.16e\n";
   printf("wcm4s: creating ascii file %s...\n",nom);
   f=fopen(nom,"w");
   if (f == NULL) {printf("wcm4s: can't open file %s\n",nom); exit(1);}
   l=0;
   for (i=ni; i<=nf; i++)
   {
      m=ntph4s(i);
      for (j=0; j<m; j++)
      {
         llex4s(j,k,i);
         a=h[i][j].re;
         b=h[i][j].im;
         c=fabs(a);
         d=fabs(b);
         if (d > tol) printf(fm1,k[0],k[1],k[2],k[3],a,b);
         if (c > tol)
            {
               e=fprintf(f,fm2,k[0],k[1],k[2],k[3],a);
               if (e == ERR)
                  {
                     fclose(f);
                     iomsg_ea2("wcm4s",l); /* this will stop the program */
                  }
               ++l;
            }
      }
   }
   printf("wcm4s: total number of (written) monomials: %d\n",l);
   printf("wcm4s: closing %s\n",nom);
   fclose(f);
   return;
}
void rrea4s(double **c, integer ni, integer nf, char *nom)
/*
this is to read a real expansion stored in the ascii file nom.
warning: before reading the file, the routine will fill with zeros
the expansion c, from degree 0 (0, not ni) to nf. if you don't
want this, you should change the marked loop below.

parameters:
c:   read series (output).
ni:  first degree to be read.
nf:  last degree to be read.
nom: name of the ascii file containing the expansion (input).
*/
{
   double a;
   integer i,j,k[6],l,m;
   char *fmt;
   FILE *f;
   fmt="%d %d %d %d %le";
   printf("rrea4s: opening ascii file %s...\n",nom);
   f=fopen(nom,"r");
   if (f == NULL) {printf("rrea4s: can't open file %s\n",nom); exit(1);}
   for (i=0; i<=nf; i++) /* this is the loop that zeroes c */
   {
      m=ntph4s(i);
      for (j=0; j<m; j++) c[i][j]=0.e0;
   }
   l=0;
   while(fscanf(f,fmt,k,k+1,k+2,k+3,&a) != EOF)
   {
      i=k[0]+k[1]+k[2]+k[3];
      k[4]=0; k[5]=0;
      if ((i >= ni) && (i <= nf))
         {
            j=exll6s(k,i);
            c[i][j]=a;
            ++l;
         }
   }
   printf("rrea4s: total number of (read) monomials: %d\n",l);
   printf("rrea4s: closing %s\n",nom);
   fclose(f);
   return;
}
void rreb4s(double **c, integer ni, integer nf, char *nom)
/*
this is to read a real expansion stored in the binary file nom.

parameters:
c:   read series (output).
ni:  first degree to be read.
nf:  last degree to be read.
nom: name of the binary file containing the expansion (input).
*/
{
   integer i,j,k,l,m,inf[4];
   double w;
   FILE *f;
   printf("rreb4s: opening binary file %s...\n",nom);
   f=fopen(nom,"rb");
   if (f == NULL) {printf("rreb4s: can't open file %s\n",nom); exit(1);}
   j=fread(inf,sizeof(integer),4,f); /* this is the head of the file */
   if (j < 4) {puts("rreb4s error 1."); exit(1);}
   if ((inf[0] != 4) || (inf[1] != 2)) iomsg_eb02("rreb4s",nom,inf[0],inf[1]);
   iomsg_eb01("rreb4s",nom,ni,nf,inf[2],inf[3]);
/*
   the procedure used to read the file is the following: first we set
   the array c (for degrees between ni and nf) to zero and then we read
   the monomials of the file whose degree belong to the interval [ni,nf].
   in this way, all the monomials that are not in the file are considered
   to be zero.
*/
   for (i=ni; i<=nf; i++)
   {
      m=ntph4s(i);
      for (j=0; j<m; j++) c[i][j]=0;
   }
   l=0;
   for (i=inf[2]; i<=inf[3]; i++)
   {
      if (i > nf) break;
      if (i >= ni)
         {
            m=ntph4s(i);
            j=fread(c[i],sizeof(double),m,f);
            if (j < m) {puts("rreb4s: reading error 1"); exit(1);}
            l += j;
         }
         else
         {
            m=ntph4s(i);
            for (j=0; j<m; j++)
            {
               k=fread(&w,sizeof(double),m,f);
               if (k < 1) {puts("rreb4s: reading error 2"); exit(1);}
            }
         }
   }
   printf("rreb4s: total number of (read) monomials: %d\n",l);
   printf("rreb4s: closing %s\n",nom);
   fclose(f);
   return;
}
