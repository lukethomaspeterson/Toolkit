/*
this file contains routines to evaluate real power expansions, of 6
variables and with the simmetries 6s and 6p. they are used to perform
changes of variables.
*/

#include <stdio.h>
#include <stdlib.h>

#include "arit-r.h"
#include "mp4s.h"
#include "mp6s.h"
#include "mp6p.h"

static renum *v[6];
static char flag='n';

void ini_erpe6(renum x[6], integer n)
/*
this is to initialize the table v, to be used by routines erpe6s
and/or erpe6p. the table will contain, after the initialization, the
different powers of each of the components of x. once the point x has
been transformed by the change the user must call end_erpe6 to free
table v. to transform a new point a new call to to ini_erpe6 must be
done to initialize v again. if end_free has not been called, the next
call to ini_erpe6 will not work.

parameters:
x: point to be transformed by the change (input).
n: degree up to which the powers of each x[i] will be computed.
*/
{
   integer i,j;
/*
   let us check that v is not initialized
*/
   if (flag != 'n')
      {
         puts("ini_erpe6s: call end_erpe6 before calling ini_erpe6s again");
         exit(1);
      }
   flag='y';
/*
   we allocate room for v
*/
   for (i=0; i<6; i++)
   {
      v[i]=(renum*)malloc((n+1)*sizeof(renum));
      if (v[i] == NULL) {puts("ini_erpe6: out of memory"); exit(1);}
   }
/*
   v is filled
*/
   for (i=0; i<6; i++)
   {
      v[i][0]=1;
      for (j=1; j<=n; j++) v[i][j]=x[i]*v[i][j-1];
   }
   return;
}
void end_erpe6(void)
/*
this is to free the table v allocated by ini_erpe6.
*/
{
   integer i;
   flag='n';
   for (i=0; i<6; i++) free(v[i]);
   return;
}

void ini_erpe4(renum x[4], integer n)
/*
this is to initialize the table v, to be used by routines erpe6s
and/or erpe6p. the table will contain, after the initialization, the
different powers of each of the components of x. once the point x has
been transformed by the change the user must call end_erpe6 to free
table v. to transform a new point a new call to to ini_erpe6 must be
done to initialize v again. if end_free has not been called, the next
call to ini_erpe6 will not work.

parameters:
x: point to be transformed by the change (input).
n: degree up to which the powers of each x[i] will be computed.
*/
{
   integer i,j;
/*
   let us check that v is not initialized
*/
   if (flag != 'n')
      {
         puts("ini_erpe6s: call end_erpe6 before calling ini_erpe4s again");
         exit(1);
      }
   flag='y';
/*
   we allocate room for v
*/
   for (i=0; i<4; i++)
   {
      v[i]=(renum*)malloc((n+1)*sizeof(renum));
      if (v[i] == NULL) {puts("ini_erpe4: out of memory"); exit(1);}
   }
/*
   v is filled
*/
   for (i=0; i<4; i++)
   {
      v[i][0]=1;
      for (j=1; j<=n; j++) v[i][j]=x[i]*v[i][j-1];
   }
   return;
}
void end_erpe4(void)
/*
this is to free the table v allocated by ini_erpe6.
*/
{
   integer i;
   flag='n';
   for (i=0; i<4; i++) free(v[i]);
   return;
}
renum erpe6s(renum **p, integer n, renum *err)
/*
this is to evaluate a real power series of the type 6s, assuming that
it contains monomials of degrees between 1 and n. it is assumed that v
has been previously initializated by ini_erpe6, with the point we want
to transform.

parameters:
p:   series to be evaluated (input).
n:   degree up to which the series is going to be evaluated. of course,
     it is assumed that the expansion p is, at least, of degree n.
err: this is the result of evaluating the monomials of degree
     exactly n. it can be seen as an heuristic estimate of the error
     (of the truncated series p with respect to the true (and infinite)
     expansion).
returned value: the result of the evaluation.
*/
{
   renum s,u;
   integer i,j,k[6],m;
   m=ntph6s(n);
   s=0;
   for (j=0; j<m; j++)
   {
      llex6s(j,k,n);
      u=v[0][k[0]]*v[1][k[1]]*v[2][k[2]]*v[3][k[3]]*v[4][k[4]]*v[5][k[5]];
      s += p[n][j]*u;
   }
   *err=s;
   for (i=n-1; i>0; i--)
   {
      m=ntph6s(i);
      for (j=0; j<m; j++)
      {
         llex6s(j,k,i);
         u=v[0][k[0]]*v[1][k[1]]*v[2][k[2]]*v[3][k[3]]*v[4][k[4]]*v[5][k[5]];
         s += p[i][j]*u;
      }
   }
   return(s);
}
renum erpe4s(renum **p, integer n, renum *err)
/*
this is to evaluate a real power series of the type 6s, assuming that
it contains monomials of degrees between 1 and n. it is assumed that v
has been previously initializated by ini_erpe6, with the point we want
to transform.

parameters:
p:   series to be evaluated (input).
n:   degree up to which the series is going to be evaluated. of course,
     it is assumed that the expansion p is, at least, of degree n.
err: this is the result of evaluating the monomials of degree
     exactly n. it can be seen as an heuristic estimate of the error
     (of the truncated series p with respect to the true (and infinite)
     expansion).
returned value: the result of the evaluation.
*/
{
   renum s,u;
   integer i,j,k[4],m;
   m=ntph4s(n);
   s=0;
   for (j=0; j<m; j++)
   {
      llex4s(j,k,n);
      u=v[0][k[0]]*v[1][k[1]]*v[2][k[2]]*v[3][k[3]];
      s += p[n][j]*u;
   }
   *err=s;
   for (i=n-1; i>0; i--)
   {
      m=ntph4s(i);
      for (j=0; j<m; j++)
      {
         llex4s(j,k,i);
         u=v[0][k[0]]*v[1][k[1]]*v[2][k[2]]*v[3][k[3]];
         s += p[i][j]*u;
      }
   }
   return(s);
}
renum erpe6p(renum **p, integer n, renum *err)
/*
this is to evaluate a real power series of the type 6p, assuming that
it contains monomials of degrees between 1 and n. it is assumed that v
has been previously initializated by ini_erpe6, with the point we want
to transform.

parameters:
p:   series to be evaluated (input).
n:   degree up to which the series is going to be evaluated. of course,
     it is assumed that the expansion p is, at least, of degree n.
err: this is the result of evaluating the monomials of degree
     exactly n. it can be seen as an heuristic estimate of the error
     (of the truncated series p with respect to the real and infinite
     expansion).
returned value: the result of the evaluation.
*/
{
   renum s,u;
   integer i,j,k[6],m;
   m=ntph6p(n);
   s=0;
   for (j=0; j<m; j++)
   {
      llex6p(j,k,n);
      u=v[0][k[0]]*v[1][k[1]]*v[2][k[2]]*v[3][k[3]]*v[4][k[4]]*v[5][k[5]];
      s += p[n][j]*u;
   }
   *err=s;
   for (i=n-1; i>0; i--)
   {
      m=ntph6p(i);
      for (j=0; j<m; j++)
      {
         llex6p(j,k,i);
         u=v[0][k[0]]*v[1][k[1]]*v[2][k[2]]*v[3][k[3]]*v[4][k[4]]*v[5][k[5]];
         s += p[i][j]*u;
      }
   }
   return(s);
}
