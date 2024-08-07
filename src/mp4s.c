/*
this file contais the basic routines to handle polynomials
of 4 variables such that the sum of the exponents of the last
two variables is always even.
*/

#include <stdio.h>
#include <stdlib.h>

#include "msgs.h"

static integer nor=0;
static integer *numo;
static natural **clmo;
/*
internal variables:
   nor:  maximum degree allowed.
   numo: numo[i] contains the number of monomials of degree i
   clmo: clmo[i][j] contains (encoded) the multiindex corresponding
         to the monomial number j inside the monomials of degree i.
*/
integer imp4s(integer nr)
/*
this routine initializes the tables used by the manipulator. it has to
be called before using the manipulator.

parameters:
nr: maximum degree we are going to work with.
    it can not be greater than 255.

returned value: number of kbytes allocated by the internal tables.
*/
{
   void prxk4s(natural k[]);
   integer i,j,nm,nt;
   unsigned long int mem;
   natural k[4];
   if (sizeof(integer) < 4) {puts("imp4s: integers need 4 bytes"); exit(1);}
   if (nr > 255) {puts("imp4s: maximum degree allowed is 255"); exit(1);}
   if (nor != 0) /* this means that imp4s is already initialized !!! */
      {
         if (nr <= nor)
            {
               mpmsg_i1("imp4s",nor);
               return(0);
            }
         mpmsg_i2("imp4s",nor,nr,"amp4s"); /* this will stop the program */
      }
   nor=nr;
   mem=(nr+1)*sizeof(integer); /* this is to count the bytes allocated */
   numo=(integer*)malloc((nr+1)*sizeof(integer));
   if (numo == NULL) {puts("imp4s error. no memory (1)."); exit(1);}
   for (i=0; i<=nr; i++)
   {
      nt=0;
      for (j=i; j>=0; j-=2) nt += (j+1)*(i-j+1);
      numo[i]=nt;
   }
   mem += (nr+1)*sizeof(natural*);
   clmo=(natural**)malloc((nr+1)*sizeof(natural*));
   if (clmo == NULL) {puts("imp4s error. no memory (2)."); exit(1);}
   for (i=0; i<=nr; i++)
   {
      nm=numo[i];
      mem += nm*sizeof(natural);
      clmo[i]=(natural*)malloc(nm*sizeof(natural));
      if (clmo[i] == NULL) {puts("imp4s error. no memory (3)."); exit(1);}
   }
   for (i=0; i<=nr; i++)
   {
      k[0]=i; k[1]=k[2]=k[3]=0;
      clmo[i][0]=i;
      nm=numo[i];
      for (j=1; j<nm; j++)
      {
         prxk4s(k);
         clmo[i][j]=k[0]+256L*k[1]+65536L*k[2]+16777216L*k[3];
      }
   }
   mem /= 1024;
   return(mem);
}
void amp4s(void)
/*
this routine frees the memory allocated by imp4s. it should be called
after using the manipulator mp4s.
*/
{
   integer i;
   if (nor == 0)
      {
         mpmsg_a1("amp4s");
         return;
      }
   for (i=0; i<=nor; i++) free(clmo[i]);
   free(clmo);
   free(numo);
   nor=0;
   return;
}
void llex4s(integer lloc, integer k[], integer no)
/*
   this routine computes the exponent of the monomial that is in
   the place number lloc of a given degree.

   parameters:
   lloc: relative place (with respect to degree no) we are interested
         in (input).
   k:    exponent of the relative place lloc (output).
   no:   order we are actually working with (input).
*/
{
   natural n;
   if (lloc > numo[no]) {puts("llex4s error."); exit(1);}
   n=clmo[no][lloc];
   k[0]=n%256;
   n/=256;
   k[1]=n%256;
   n/=256;
   k[2]=n%256;
   k[3]=n/256;
   return;
}
integer exll4s(integer k[],integer no)
/*
   this routine computes the place corresponding to the multiindex k.

   parameters:
   k:    multiindex (input).
   no:   degree we are actually working with (input).

   returned value: place corresponding to k.
*/
{
   integer i,m,n,lloc,n2;
   n=k[0]+k[1]+k[2]+k[3];
   if (n != no) {puts("exll4s error 1."); exit(1);}
   m=k[3];
   //if (m%2 != 0) {puts("exll4s error 2."); exit(1);}
   n2=n-m;
   lloc=k[3]*(n2+1);
   for (i=0; i<m; i+=2) lloc += (i+1)*(n-i+1);
   lloc += k[1];
   return(lloc);
}
integer ntph4s(integer no)
/*
   this routine returns the number of monomials of degree no.

   parameters:
   no: order we are interested in (input).

   returned value: number of monomials of order no.
*/
{
   if (no > nor) {puts("ntph4s: error."); exit(1);}
   return(numo[no]);
}
void prxk4s(natural k[])
/*
   given a multiindex k, this routine computes the next one
   according the lexicographic order. We recall that k[2]+k[3]
   must be even.
   warning: if the last multiindex of a given order is used,
   this routine could not work properly.

   parameters:
   k: array of 4 components containing the multiindex. it is
      overwritten on exit (input and output).
*/
{
   if (k[0] != 0) {k[0]--; k[1]++; return;}
   if (k[2] != 0) {k[2]--; k[3]++; k[0]=k[1]; k[1]=0; return;}
   if (k[1] <= 1) {puts("prxk4s: error 1."); exit(1);}
   k[2]=k[3]+2;
   k[3]=0;
   k[0]=k[1]-2;
   k[1]=0;
   return;
}
