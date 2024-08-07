#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "arit-r.h"

typedef int integer;

void litr5(renum a[6][6], renum mu, renum x[6], renum y[6], integer dir)
/*
this is to transform points by means of the linear change that puts
the hamiltonian in linear (and real) normal form around l5. the change
also includes the translation that moves the origin from the centre of
masses (in the synodical system) to the triangular point l5. you can
send points from normal form coordinates to synodical ones (dir=1), or
from synodical ones to the ones of the linear normal form (dir=-1).

parameters:
a:   change used to diagonalize h2 (input).
x:   initial point (input).
y:   final point (output).
dir: direction of the transformation. 1 means that the change is given
     by matrix a, -1 means that it is given by the inverse of a. as a
     is a symplectic matrix, we will compute its inverse using this
     property.
*/
{
   static integer k=0;
   static renum r;
   integer i;
   renum c,v[6];
   if (k == 0)
      {
         k=1;
         r=sqrt(3.e0)/2.e0;
      }
   c=mu-0.5;
   switch (dir)
   {
      case 1:
/*
         we first multiply by a
*/
         for (i=0; i<6; i++)
            y[i]=a[i][0]*x[0]+a[i][1]*x[1]+a[i][2]*x[2]
                +a[i][3]*x[3]+a[i][4]*x[4]+a[i][5]*x[5];
/*
         and then we translate the origin to the centre of masses
*/
         y[0] += c;
         y[1] -= r;
         y[2] += r;
         y[3] += c;
         break;
      case -1:
/*
         first, we translate the origin to l5
*/
         y[0]=x[0]-c;
         y[1]=x[1]+r;
         y[2]=x[2]-r;
         y[3]=x[3]-c;
         y[4]=x[4];
         y[5]=x[5];
/*
         next we multiply by the inverse of a. as a is a symplectic
         matrix, we do this in four steps:
         1) multiplication by the matrix of the symplectic form
*/
         v[0]=y[1]; v[1]=-y[0];
         v[2]=y[3]; v[3]=-y[2];
         v[4]=y[5]; v[5]=-y[4];
/*
         2) product by the transpose of a
*/
         for (i=0; i<6; i++)
            y[i]=a[0][i]*v[0]+a[1][i]*v[1]+a[2][i]*v[2]
                +a[3][i]*v[3]+a[4][i]*v[4]+a[5][i]*v[5];
/*
         3) product by the matrix of the symplectic form
*/
         v[0]=y[1]; v[1]=-y[0];
         v[2]=y[3]; v[3]=-y[2];
         v[4]=y[5]; v[5]=-y[4];
/*
         4) change of sign
*/
         for (i=0; i<6; i++) y[i]=-v[i];
         break;
      default:
         printf("litr5 error: dir must be 1 or -1, and it is %d\n",dir);
         exit(1);
   }
   return;
}
