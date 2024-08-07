/*
here we define the file names to be used along the program. some of
these files are ascii files, some others are binary. see the
documentation for more details.
*/

#define CTL "nf.ctl" /* to write the parameters used */
#define CVL "nf.cvl" /* to write the linear change that diagonalizes h2 */
#define GEN "nf.gen" /* to write the generating function */
#define CMP "nf.cmp" /* complex normal form (optional) */
#define RES "nf.res" /* real normal form */
#define WRK "nf.tmp" /* temporary file. erased on successful exit. */
#define CD1 "cvnf.1"
#define CD2 "cvnf.2"
#define CD3 "cvnf.3" /*  changes of variables (from normal form  */
#define CD4 "cvnf.4" /*  to initial coordinates)                 */
#define CD5 "cvnf.5"
#define CD6 "cvnf.6"
#define CI1 "cvnfi.1"
#define CI2 "cvnfi.2"
#define CI3 "cvnfi.3" /*  changes of variables (from initial to  */
#define CI4 "cvnfi.4" /*  normal form coordinates)               */
#define CI5 "cvnfi.5"
#define CI6 "cvnfi.6"
#define HEX "Hseries.dat" /* to write the series expansion of Hamiltonian */
