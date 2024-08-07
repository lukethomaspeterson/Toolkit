/*
this is to print the error messages coming from mp??.c files.
*/

#include <stdio.h>
#include <stdlib.h>

#if !defined(INTEGER)
#define INTEGER

typedef int integer;

#endif


void mpmsg_i1(char *who, integer nor)
{
   puts("**************************************************");
   printf("%s warning message:\n",who);
   printf("%s is already initialized to degree %d\n",who,nor);
   printf("now you want to initialize it again to the same degree.\n");
   printf("action taken: this call to %s is simply ignored.\n",who);
   puts("**************************************************");
   return;
}
void mpmsg_i2(char *who, integer nor, integer nr, char *amp)
{
   printf("%s error message:\n",who);
   printf("%s is already initialized to degree %d\n",who,nor);
   printf("and now you want to initialize it to degree %d.\n",nr);
   printf("you must call routine %s first.\n",amp);
   printf("action taken: program aborted\n");
   exit(1);
}
void mpmsg_a1(char *who)
{
   puts("**************************************************");
   printf("%s warning message:\n",who);
   printf("no memory to free\n");
   printf("action taken: this call is simply ignored.\n");
   puts("**************************************************");
   return;
}
void iomsg_ea2(char *who, integer n)
{
   printf("%s error message:\n",who);
   printf("writing error (disk full?)\n");
   printf("number of lines written: %d\n",n);
   printf("action taken: program aborted\n");
   exit(1);
}
void iomsg_eb01(char *who, char *file, integer ni, integer nf,
               integer nif, integer nff)
{
   if ((ni >= nif) && (nf <= nff)) return;
   puts("**************************************************");
   printf("%s warning message:\n",who);
   printf("you want to read from degree %d to degree %d\n",ni,nf);
   printf("the file contains from degree %d to degree %d\n",nif,nff);
   printf("action taken: the missing monomials are considered zero.\n");
   puts("**************************************************");
   return;
}
void iomsg_eb02(char *who, char *file, integer nvf, integer s)
{
   printf("%s error message:\n",who);
   printf("file %s contains an expansion of %d variables\n",file,nvf);
   printf("with the simmetry code %d\n",s);
   printf("action taken: program aborted\n");
   exit(1);
}
void iomsg_eb2(char *who, char *file)
{
   printf("%s error message:\n",who);
   printf("writing error (disk full?)\n");
   printf("file: %s\n",file);
   printf("action taken: program aborted\n");
   exit(1);
}
