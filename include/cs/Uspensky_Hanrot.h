/*
   Uspensky's algorithm : isolation of real roots of a polynomial
   based on Descartes' rule, following the procedure described by
   Rouillier & Zimmermann

   Requires GMP.

   written by G. Hanrot (04/99, 10/02), with some improvements by
   F. Rouillier.

   Thanks to S. Lazard and S. Petitjean for feedback and fixes.

   Tested on the following platforms :
   AMD K7, GMP 4.1, LINUX 2.4.19, gcc 2.95.3
   Alpha ev6, GMP 4.1, Tru64 v4.0, gcc 2.95.3
*/
#ifndef USPENSKY_HANROT_H
#define USPENSKY_HANROT_H

#include <gmp.h>
#include <stdio.h>

/*
   Recommended in a "general purpose context" if implicit bounds
   for intervals are allowed.
*/
//#define POWER_HACK

typedef struct
{
  mpz_t c;
  long k;
  unsigned int isexact;
#ifdef POWER_HACK
  unsigned int sign; /* For roots of the type -sqrt(x) \neq sqrt(-x)... */
#endif
} interval;

#ifdef __cplusplus
extern "C" {
#endif

interval *
Uspensky(mpz_t *Q, unsigned long deg, unsigned int *nbroot);

void
affiche_root(FILE *stream, interval z);

#ifdef __cplusplus
}
#endif

#endif // USPENSKY_HANROT_H
