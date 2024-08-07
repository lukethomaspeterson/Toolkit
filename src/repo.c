/*
this is to read the ascii files with the points to be transformed by a
change of coordinates. the format of the file must be:

n d
x1 x2 x3 x4 x5 x6
y1 y2 y3 y4 y5 y6
.
.
.

or

n d
x1 x2 x3 x4
y1 y2 y3 y4
.
.
.

where n is the number of points to be read from the file, and d is the
number of coordinates for each point (it must be 6 for the first
example and 4 for the second one).
*/

#include <stdio.h>
#include <stdlib.h>

#include "arit-r.h"

renum** repo(char *nom, int d, int *np)
/*
this is to read the file with the points we are going to transform.
we will assume that the points we are going to read are double
constants, and that renum is an alias of double.

parameters:
nom: name of the file with the points (input).
d:   dimension of the points, it must be 4 or 6 (input).
np:  number of read points (output).

returned value: table filled with the read points. the first index
     corresponds to the points, and the second one to the components
     of each point.
*/
{
   double **p;
   int i,j,n,nd;
   FILE *f;
   char *fmt4,*fmt6;
   if ((d != 4) && (d != 6)) {printf("repo error 1: d=%d\n",d); exit(1);}
   fmt4="%le %le %le %le";         /* format to read the points */
   fmt6="%le %le %le %le %le %le"; /* format to read the points */
   f=fopen(nom,"r");
   if (f == NULL) {printf("repo: cannot open %s\n",nom); exit(1);}
   fscanf(f,"%d %d",&n,&nd);
   if (nd != d) {printf("repo error 2: d=%d nd=%d\n",d,nd); exit(1);}
   p=(renum**)malloc(n*sizeof(renum*));
   if (p == NULL) {puts("repo: out of memory (1)"); exit(1);}
   for (i=0; i<n; i++)
   {
      p[i]=(renum*)malloc(d*sizeof(renum));
      if (p[i] == NULL) {puts("repo: out of memory (2)"); exit(1);}
   }
   if (d == 4)
      {
         for (j=0; j<n; j++)
            fscanf(f,fmt4,p[j],p[j]+1,p[j]+2,p[j]+3);
      }
      else
      {
         for (j=0; j<n; j++)
            fscanf(f,fmt6,p[j],p[j]+1,p[j]+2,p[j]+3,p[j]+4,p[j]+5);
      }
   fclose(f);
   *np=n;
   return(p);
}
