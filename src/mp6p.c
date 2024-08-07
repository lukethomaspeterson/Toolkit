/*
this file contais the basic routines to handle polynomials
of 6 variables such that the sum of the exponents of the last
two variables is always odd.
*/

#include <stdio.h>
#include <stdlib.h>

#include "msgs.h"

static integer nor=0;
static natural **clmo;
static integer **psi;
static integer *numo;
/*
internal variables:
   nor:  maximum degree allowed.
   clmo: clmo[i][j] contains (encoded) the multiindex corresponding
         to the monomial number j inside the monomials of degree i.
   psi:  psi[i][j] contains $\psi_{i}(j)$.
   numo: numo[i] contains the number of monomials of degree i.

remark: psi[6][i] is not numo[i], because of the symmetry (psi[6][i]
   is the total number of monomials without taking into account the
   symmetry, and numo[i] takes it into account). in fact, psi[i][j] is
   only computed for i=2,3,4.
*/

integer imp6p(integer nr)
/*
this routine initializes the tables used by the manipulator. it has to
be called before using the manipulator.

parameters:
nr: maximum degree we are going to work with.
    it can not be greater than 63.

returned value: number of kbytes allocated by the internal tables.
*/
{
   void prxk6p(natural k[]);
   integer i,j,nm,l;
   unsigned long int mem;
   natural k[6];
   if (sizeof(integer) < 4) {puts("imp6p: integers need 4 bytes"); exit(1);}
   if (nr > 63) {puts("imp6p: maximum degree allowed is 63"); exit(1);}
   if (nor != 0) /* this means that imp6p is already initialized !!! */
      {
         if (nr <= nor)
            {
               mpmsg_i1("imp6p",nor);
               return(0);
            }
         mpmsg_i2("imp6p",nor,nr,"amp6p"); /* this will stop the program */
      }
   nor=nr;
   mem=3*sizeof(integer); /* this is to count the bytes allocated */
   psi=(integer**)malloc(3*sizeof(integer*));
   if (psi == NULL) {puts("imp6p error. no memo (1)."); exit(1);}
   psi -= 2;
   for (i=2; i<=4; i++)
   {
      psi[i]=(integer*)malloc((nr+1)*sizeof(integer));
      if (psi[i] == NULL) {puts("imp6p error. no memo (2)."); exit(1);}
      mem += (nr+1)*sizeof(integer);
   }
   for (j=0; j<=nr; j++) psi[2][j]=j+1;
   for (i=3; i<=4; i++)
   {
      for (j=0; j<=nr; j++)
      {
         psi[i][j]=0;
         for (l=0; l<=j; l++) psi[i][j] += psi[i-1][l];
      }
   }
   mem += (nr+1)*sizeof(integer);
   numo=(integer*)malloc((nr+1)*sizeof(integer));
   if (numo == NULL) {puts("imp6p error. no memo (3)."); exit(1);}
   mem += (nr+1)*sizeof(natural*);
   clmo=(natural**)malloc((nr+1)*sizeof(natural*));
   if (clmo == NULL) {puts("imp6p error. no memo (4)."); exit(1);}
   for (i=0; i<=nr; i++)
   {
      nm=0;
      for (j=1; j<=i; j+=2) nm += psi[4][i-j]*(j+1);
      numo[i]=nm;
      mem += nm*sizeof(natural);
      clmo[i]=(natural*)malloc(nm*sizeof(natural));
      if (clmo[i] == NULL) {puts("imp6p error. no memo (5)."); exit(1);}
   }
   for (i=1; i<=nr; i++)
   {
      k[0]=i-1; k[1]=k[2]=k[3]=k[5]=0; k[4]=1;
      clmo[i][0]=262144L*k[4];
      nm=numo[i];
      for (j=1; j<nm; j++)
      {
         prxk6p(k);
         clmo[i][j]=k[1]+64L*k[2]+4096L*k[3]+262144L*k[4]+16777216L*k[5];
      }
   }
   mem /= 1024;
   return(mem);
}
void amp6p(void)
/*
this routine frees the memory allocated by imp6p. it should be called
after using the manipulator mp6p.
*/
{
   integer i;
   if (nor == 0)
      {
         mpmsg_a1("amp6p");
         return;
      }
   for (i=0; i<=nor; i++) free(clmo[i]);
   free(clmo);
   free(numo);
   for (i=2; i<=4; i++) free(psi[i]);
   psi += 2;
   free(psi);
   nor=0;
   return;
}
void llex6p(integer lloc, integer k[], integer no)
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
   integer m;
   if (lloc > numo[no]) {puts("llex6p error."); exit(1);}
   n=clmo[no][lloc];
   k[1]=n%64;
   m=k[1];
   n/=64;
   k[2]=n%64;
   m+=k[2];
   n/=64;
   k[3]=n%64;
   m+=k[3];
   n/=64;
   k[4]=n%64;
   m+=k[4];
   k[5]=n/64;
   m+=k[5];
   k[0]=no-m;
   return;
}
integer exll6p(integer k[],integer no)
/*
   this routine computes the place corresponding to the multiindex k.

   parameters:
   k:    multiindex (input).
   no:   degree we are actually working with (input).

   returned value: place corresponding to k.
*/
{
   integer i,n,m,lloc,n4,n3,n2;
   n=k[0]+k[1]+k[2]+k[3]+k[4]+k[5];
   if (n != no) {puts("exll6p error 1."); exit(1);}
   m=k[4]+k[5];
   if (m%2 != 1) {puts("exll6p error 2."); exit(1);}
   n4=n-m;
   lloc=k[5]*psi[4][n4];
   for (i=1; i<m; i+=2) lloc += (i+1)*psi[4][n-i];
   n3=n4-k[3];
   for (i=n3+1; i<=n4; i++) lloc += psi[3][i];
   n2=n3-k[2];
   for (i=n2+1; i<=n3; i++) lloc += psi[2][i];
   lloc += k[1];
   return(lloc);
}
integer ntph6p(integer no)
/*
   this routine returns the number of monomials of degree no.

   parameters:
   no: order we are interested in (input).

   returned value: number of monomials of order no.
*/
{
   if (no > nor) {puts("ntph6p: error."); exit(1);}
   return(numo[no]);
}
void prxk6p(natural k[])
/*
   given a multiindex k, this routine computes the next one
   according to a product-lexicographic order. We recall that
   k[4]+k[5] must be odd.
   warning: if the last multiindex of a given order is used,
   this routine could not work properly.

   parameters:
   k: array of 6 components containing the multiindex. it is
      overwritten on exit (input and output).
*/
{
   if (k[0] != 0) {k[0]--; k[1]++; return;}
   if (k[1] != 0) {k[0]=k[1]-1; k[1]=0; k[2]++; return;}
   if (k[2] != 0) {k[0]=k[2]-1; k[2]=0; k[3]++; return;}
   if (k[4] != 0) {k[4]--; k[5]++; k[0]=k[3]; k[3]=0; return;}
   if (k[3] <= 1) {puts("prxk6p: error 1."); exit(1);}
   k[4]=k[5]+2;
   k[5]=0;
   k[0]=k[3]-2;
   k[3]=0;
   return;
}
