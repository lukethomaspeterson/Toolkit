#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* uneix(char *a, char *b)
/*
this is to merge two strings.

parameters:
a: first string (input).
b: second string (input).
returned value: the result of putting string b at the end of string a
*/
{
   char *c;
   int n;
   n=strlen(a)+strlen(b)+1;
   c=(char*)malloc(n*sizeof(char));
   if (c == NULL) {puts("uneix error: out of memory"); exit(1);}
   c=strcpy(c,a);
   c=strcat(c,b);
   return(c);
}
