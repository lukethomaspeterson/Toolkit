/*
this file contains a routine that writes the linear change of
variables that put h2 in real normal form. it is assumed that the
standard double precision arithmetic is used.

here you will also find two small routines that read and write a small
ascii file with the parameters used in the actual run of the program.
*/

#include <stdio.h>
#include <stdlib.h>

typedef int integer;
typedef double renum;

void wcvl(renum c[6][6], char *nom)
/*
this is to write (in an ascii file) the real linear change that
puts h2 in real normal form. it is needed to send points from the normal
form coordinates to the initial ones.

parameters:
c:   matrix with the linear change (input).
nom: name of the file where the change will be stored (input).
*/
{
   integer i,j;
   FILE *f;
   f=fopen(nom,"w");
   if (f == NULL) {printf("wcvl: cannot open file %s\n",nom); exit(1);}
   for (i=0; i<6; i++)
   {
      for (j=0; j<6; j++)
         fprintf(f,"%2d %2d %24.16e\n",i,j,c[i][j]);
   }
   fclose(f);
   return;
}
void rcvl(renum c[6][6], char *nom)
/*
this is to read (from an ascii file previously written by wcvl) the
real linear change that puts h2 in real normal form.

parameters:
c:   matrix with the linear change (output).
nom: name of the file where the change is stored (input).
*/
{
   integer i,j,k1,k2;
   renum z;
   FILE *f;
   f=fopen(nom,"r");
   if (f == NULL) {printf("rcvl: cannot open file %s\n",nom); exit(1);}
   for (i=0; i<6; i++)
   {
      for (j=0; j<6; j++)
      {
         fscanf(f,"%d %d %le",&k1,&k2,&z);
         if (i != k1) puts("rcvl: warning (index mismatch) 1.");
         if (j != k2) puts("rcvl: warning (index mismatch) 2.");
         c[i][j]=z;
      }
   }
   fclose(f);
   return;
}
void wctl123(char *ctl, integer li, integer n, renum mu, renum gam)
/*
this routine writes, in an ascii file, the parameters used in the
actual run.

parameters:
ctl: name of the file where the parameters have to be stored (input).
li:  libration point.
n:   degree of the expansions.
mu:  mass parameter of the rtbp.
gam: solution of the corresponding euler quintic equation.
*/
{
   FILE *f;
   f=fopen(ctl,"w");
   if (f == NULL) {printf("wctl123: cannot open file %s\n",ctl); exit(1);}
   fprintf(f,"%4d\n%4d\n%24.16e\n%24.16e\n",li,n,mu,gam);
   fclose(f);
   return;
}
void rctl123(char *ctl, integer *li, integer *n, renum *mu, renum *gam)
/*
this routine reads, from an ascii file (previously written by wctl),
the parameters used in the actual run.

parameters:
ctl: name of the file containing the parameters (input).
li:  libration point (output).
n:   degree of the expansions (output).
mu:  mass parameter of the rtbp (output).
gam: solution of the corresponding euler quintic equation (output).
*/
{
   FILE *f;
   f=fopen(ctl,"r");
   if (f == NULL) {printf("rctl123: cannot open file %s\n",ctl); exit(1);}
   fscanf(f,"%d",li);
   fscanf(f,"%d",n);
   fscanf(f,"%le",mu);
   fscanf(f,"%le",gam);
   fclose(f);
   return;
}
void wctl5(char *ctl, integer n, renum mu)
/*
this routine writes, in an ascii file, the parameters used in the
actual run (l5 case).

parameters:
ctl: name of the file where the parameters have to be stored (input).
n:   degree of the expansions.
mu:  mass parameter of the rtbp.
*/
{
   FILE *f;
   f=fopen(ctl,"w");
   if (f == NULL) {printf("wctl5: cannot open file %s\n",ctl); exit(1);}
   fprintf(f,"%4d\n%24.16e\n",n,mu);
   fclose(f);
   return;
}
void rctl5(char *ctl, integer *n, renum *mu)
/*
this routine reads, from an ascii file, the parameters used in the
actual run (l5 case).

parameters:
ctl: name of the file where the parameters have to be found (input).
n:   degree of the expansions (output).
mu:  mass parameter of the rtbp (output).
*/
{
   FILE *f;
   f=fopen(ctl,"r");
   if (f == NULL) {printf("rctl5: cannot open file %s\n",ctl); exit(1);}
   fscanf(f,"%d",n);
   fscanf(f,"%le",mu);
   fclose(f);
   return;
}

