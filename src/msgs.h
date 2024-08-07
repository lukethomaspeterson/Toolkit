#if !defined(INTEGER)
#define INTEGER

typedef int integer;
typedef unsigned int natural;

#endif

void mpmsg_i1(char *who, integer nor);
void mpmsg_i2(char *who, integer nor, integer nr, char *amp);
void mpmsg_a1(char *who);

void iomsg_ea2(char *who, integer n);
void iomsg_eb01(char *who, char *file, integer ni, integer nf,
               integer nif, integer nff);
void iomsg_eb02(char *who, char *file, integer nvf, integer s);
void iomsg_eb2(char *who, char *file);
