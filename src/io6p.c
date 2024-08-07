/*
this file contains high-level i/o functions (in C) for the polynomials
handled by mp6p. note that the definition of 'complex' is done 'by
hand', instead of including the file arit-c.h. so, a complex is a
structure containing two doubles.
*/

#include <stdio.h>
#include <stdlib.h>

#include "mp6p.h"
#include "msgs.h"

typedef struct {
   double re;
   double im;
} complex;

#define ERR -1 /* value returned by fprintf on failure */

void wea6p(complex **h, integer ni, integer nf, char *nom)
/*
this is to write, in ascii format, the expansion h in the file nom.
the file is also created by this routine. if the file already exists,
the previous one is deleted.

parameters:
h:   expansion to be written (input).
ni:  degree to start writting the hamiltonian.
nf:  last degree to be written.
nom: string with the name of the file to write h (input).
*/
{
   integer i,j,k[6],l,m,c;
   FILE *f;
   printf("wea6p: creating ascii file %s...\n",nom);
   f=fopen(nom,"w");
   if (f == NULL) {printf("wea6p: can't open file %s\n",nom); exit(1);}
   l=0;
   for (i=ni; i<=nf; i++)
   {
      m=ntph6p(i);
      for (j=0; j<m; j++)
      {
         llex6p(j,k,i);
         c=fprintf(f,"%2d %2d %2d %2d %2d %2d %25.16e %25.16e\n",
             k[0],k[1],k[2],k[3],k[4],k[5],h[i][j].re,h[i][j].im);
         if (c == ERR)
            {
               fclose(f);
               iomsg_ea2("wea6p",l); /* this will stop the program */
            }
         ++l;
      }
   }
   printf("wea6p: total number of (written) monomials: %d\n",l);
   printf("wea6p: closing %s\n",nom);
   fclose(f);
   return;
}
void rea6p(complex **h, integer ni, integer nf, char *nom)
/*
this is to read a real expansion stored in the ascii file nom.
warning: before reading the file, the routine will fill with zeros
the expansion h, from degree 0 (0, not ni) to nf. if you don't
want this, you should change the marked loop below.

parameters:
h:   read series (output).
ni:  first degree to be read.
nf:  last degree to be read.
nom: name of the ascii file containing the expansion (input).
*/
{
   double u,v;
   integer i,j,k[6],m;
   char *fmt;
   FILE *f;
   fmt="%d %d %d %d %d %d %le %le";
   printf("rea6p: opening ascii file %s...\n",nom);
   f=fopen(nom,"r");
   if (f == NULL) {printf("rea6p: can't open file %s\n",nom); exit(1);}
   for (i=0; i<=nf; i++) /* this is the loop that zeroes h */
   {
      m=ntph6p(i);
      for (j=0; j<m; j++) {h[i][j].re=0.e0; h[i][j].im=0.e0;}
   }
   m=0;
   while(fscanf(f,fmt,k,k+1,k+2,k+3,k+4,k+5,&u,&v) != EOF)
   {
      i=k[0]+k[1]+k[2]+k[3]+k[4]+k[5];
      if ((i >= ni) && (i <= nf))
         {
            j=exll6p(k,i);
            h[i][j].re=u;
            h[i][j].im=v;
            ++m;
         }
   }
   printf("rea6p: total number of (read) monomials: %d\n",m);
   printf("rea6p: closing %s\n",nom);
   fclose(f);
   return;
}
void web6p(complex **h, integer ni, integer nf, char *nom)
/*
this is to write, in binary format, the expansion h in the file nom.
the file is also created by this routine. if the file already exists,
the previous one is deleted.

parameters:
h:   expansion to be written (input).
ni:  degree to start writting the hamiltonian.
nf:  last degree to be written.
nom: string with the name of the file to write h (input).
*/
{
   integer i,j,l,m,inf[4];
   FILE *f;
   printf("web6p: creating binary file %s...\n",nom);
   f=fopen(nom,"wb");
   if (f == NULL) {printf("web6p: can't open file %s\n",nom); exit(1);}
   inf[0]=6;  /* number of variables */
   inf[1]=1;  /* kind of simmetry */
   inf[2]=ni; /* initial degree */
   inf[3]=nf; /* final degree */
   j=fwrite(inf,sizeof(integer),4,f);
   if (j < 4) iomsg_eb2("web6p",nom); /* this will stop the program */
   l=0;
   for (i=ni; i<=nf; i++)
   {
      m=ntph6p(i);
      j=fwrite(h[i],sizeof(complex),m,f);
      if (j < m) iomsg_eb2("web6p",nom); /* this will stop the program */
      l += j;
   }
   printf("web6p: total number of (written) monomials: %d\n",l);
   printf("web6p: closing %s\n",nom);
   fclose(f);
   return;
}
void reb6p(complex **h, integer ni, integer nf, char *nom)
/*
this is to read an expansion contained in the binary file nom.

parameters:
h:   read series (output).
ni:  first degree to be read.
nf:  last degree to be read.
nom: name of the binary file containing the expansion (input).
*/
{
   integer i,j,l,m,k,inf[4];
   complex w;
   FILE *f;
   printf("reb6p: opening binary file %s...\n",nom);
   f=fopen(nom,"rb");
   if (f == NULL) {printf("reb6p: can't open file %s\n",nom); exit(1);}
   j=fread(inf,sizeof(integer),4,f); /* this is the head of the file */
   if (j < 4) {puts("reb6p error 1."); exit(1);}
   if ((inf[0] != 6) || (inf[1] != 1)) iomsg_eb02("reb6p",nom,inf[0],inf[1]);
   iomsg_eb01("reb6p",nom,ni,nf,inf[2],inf[3]);
/*
   the procedure used to read the file is the following: first we set
   the array h (for degrees between ni and nf) to zero and then we read
   the monomials of the file whose degree belong to the interval [ni,nf].
   in this way, all the monomials that are not in the file are considered
   to be zero.
*/
   for (i=ni; i<=nf; i++)
   {
      m=ntph6p(i);
      for (j=0; j<m; j++) {h[i][j].re=0; h[i][j].im=0;}
   }
   l=0;
   for (i=inf[2]; i<=inf[3]; i++)
   {
      if (i > nf) break;
      if (i >= ni)
         {
            m=ntph6p(i);
            j=fread(h[i],sizeof(complex),m,f);
            if (j < m) {puts("reb6p: reading error 1"); exit(1);}
            l += j;
         }
         else
         {
            m=ntph6p(i);
            for (j=0; j<m; j++)
            {
               k=fread(&w,sizeof(complex),m,f);
               if (k < 1) {puts("reb6p: reading error 2"); exit(1);}
            }
         }
   }
   printf("reb6p: total number of (read) monomials: %d\n",l);
   printf("reb6p: closing %s\n",nom);
   fclose(f);
   return;
}
void wpb6p(complex *p, integer n, char *nom, char flag)
/*
this is to write an homogeneous polynomial p of degree n in the
binary file nom. if flag is 'a', the polynomial is appended to the
end of the file, otherwise the file is deleted before writting p.

parameters:
p:    homogeneous polynomial (input).
n:    degree of p.
nom:  name of the file where p has to be stored (input).
flag: see remark above.
*/
{
   integer j,m;
   char *mode;
   FILE *f;
   if (flag == 'a') {mode="ab";} else {mode="wb";}
   printf("wpb6p: opening file %s in mode %s...\n",nom,mode);
   f=fopen(nom,mode);
   if (f == NULL) {printf("wpb6p: can't open file %s\n",nom); exit(1);}
   printf("wpb6p: writing degree %d\n",n);
   m=ntph6p(n);
   j=fwrite(p,sizeof(complex),m,f);
   if (j < m) iomsg_eb2("wpb6p",nom); /* this will stop the program */
   printf("wpb6p: closing %s\n",nom);
   fclose(f);
   return;
}
void rpb6p(FILE *f, complex *p, integer n)
/*
this is to read a polynomial of degree n from the binary file pointed
by f. it is assumed that this file has already been opened by the
calling routine.
warning: this routine does not perform any check on the read data.

parameters:
f: pointer to the (already opened) file to be read (input).
p: read polynomial (output).
n: degree of the polynomial to be read and of p. warning: this is
   not checked by the routine.
*/
{
   integer j,m;
   printf("rpb6p: reading degree %d\n",n);
   m=ntph6p(n);
   j=fread(p,sizeof(complex),m,f);
   if (j < m) {puts("rpb6p: error 1."); exit(1);}
   return;
}
void rrea6p(double **h, integer ni, integer nf, char *nom)
/*
this is to read a real expansion stored in the ascii file nom.
warning: before reading the file, the routine will fill with zeros
the expansion h, from degree 0 (0, not ni) to nf. if you don't
want this, you should change the marked loop below.

parameters:
h:   read series (output).
ni:  first degree to be read.
nf:  last degree to be read.
nom: name of the ascii file containing the expansion (input).
*/
{
   double u;
   integer i,j,k[6],m;
   char *fmt;
   FILE *f;
   fmt="%d %d %d %d %d %d %le";
   printf("rea6p: opening ascii file %s...\n",nom);
   f=fopen(nom,"r");
   if (f == NULL) {printf("rea6p: can't open file %s\n",nom); exit(1);}
   for (i=0; i<=nf; i++) /* this is the loop that zeroes h */
   {
      m=ntph6p(i);
      for (j=0; j<m; j++) {h[i][j]=0.e0;}
   }
   m=0;
   while(fscanf(f,fmt,k,k+1,k+2,k+3,k+4,k+5,&u) != EOF)
   {
      //i=k[0]+k[1]+k[2]+k[3];
      i=k[0]+k[1]+k[2]+k[3]+k[4]+k[5];
      if ((i >= ni) && (i <= nf))
         {
            j=exll6p(k,i);
            h[i][j]=u;
            ++m;
         }
   }
   printf("rea6p: total number of (read) monomials: %d\n",m);
   printf("rea6p: closing %s\n",nom);
   fclose(f);
   return;
}
void rreb6p(double **c, integer ni, integer nf, char *nom)
/*
this is to read a real expansion stored in the binary file nom.

parameters:
c:   read series (output).
n:   maximum degree read.
nom: name of the binary file containing the expansion (input).
*/
{
   integer i,j,l,m,k,inf[4];
   double w;
   FILE *f;
   printf("rreb6p: opening binary file %s...\n",nom);
   f=fopen(nom,"rb");
   if (f == NULL) {printf("rreb6p: can't open file %s\n",nom); exit(1);}
   j=fread(inf,sizeof(integer),4,f); /* this is the head of the file */
   if (j < 4) {puts("rreb6p error 1."); exit(1);}
   if ((inf[0] != 6) || (inf[1] != 1)) iomsg_eb02("rreb6p",nom,inf[0],inf[1]);
   iomsg_eb01("rreb6p",nom,ni,nf,inf[2],inf[3]);
/*
   the procedure used to read the file is the following: first we set
   the array c (for degrees between ni and nf) to zero and then we read
   the monomials of the file whose degree belong to the interval [ni,nf].
   in this way, all the monomials that are not in the file are considered
   to be zero.
*/
   for (i=ni; i<=nf; i++)
   {
      m=ntph6p(i);
      for (j=0; j<m; j++) c[i][j]=0;
   }
   l=0;
   for (i=inf[2]; i<=inf[3]; i++)
   {
      if (i > nf) break;
      if (i >= ni)
         {
            m=ntph6p(i);
            j=fread(c[i],sizeof(double),m,f);
            if (j < m) {puts("rreb6p: reading error 1"); exit(1);}
            l += j;
         }
         else
         {
            m=ntph6p(i);
            for (j=0; j<m; j++)
            {
               k=fread(&w,sizeof(double),m,f);
               if (k < 1) {puts("rreb6p: reading error 2"); exit(1);}
            }
         }
   }
   printf("rreb6p: total number of (read) monomials: %d\n",l);
   printf("rreb6p: closing %s\n",nom);
   fclose(f);
   return;
}
