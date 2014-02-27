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
#include <cs/Uspensky_Hanrot.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <math.h>

#define vali(x) mpz_scan1((x), 0)
#define ilog2(a) mpz_sizeinbase(a,2)
#define TOT_POS -1

/* Probably out of memory long since on nontrivial examples, but anyway. */
#define DEPTH_MAX 1024

#if 0

/*
   Recommended in a "general purpose context" if implicit bounds
   for intervals are allowed.
*/
#define POWER_HACK

#endif

static unsigned long transl_done, node_looked, half_done, max_depth, sign;
static unsigned long bound_pos_flag, bound_neg_flag, b_pos_flag, b_neg_flag;
static int b1 = 0, b2 = 0;

#ifdef POWER_HACK
static unsigned long t_power;
#endif

#if 0

typedef struct
{
  mpz_t c;
  long k;
  unsigned int isexact;
#ifdef POWER_HACK
  unsigned int sign; /* For roots of the type -sqrt(x) \neq sqrt(-x)... */
#endif
} interval;

#endif

#ifdef DEBUG
void
output(mpz_t *P, unsigned long deg)
{
  unsigned long i;

  for (i = 0; i <= deg; i++)
    {
      mpz_out_str(stderr, 10, P[i]); fprintf(stderr, "*x^%ld+", i);
    }
  fprintf(stderr, "\n");
}
#endif

int
RemoveContent(mpz_t *P, unsigned long deg)
{
  unsigned long cont, i, z;

  i = 0; while (mpz_sgn(P[i]) == 0) i++;
  cont = vali(P[i]);

  for( ; (i <= deg) && cont; i++)
    {
      if (mpz_sgn(P[i]) != 0)
	{
	  z = vali(P[i]);
	  if (z < cont) cont = z;
	}
    }

#ifdef DEBUG
  fprintf(stderr, "Removing a content of 2^%ld\n", cont);
#endif
  if (cont == 0) return 0;

  for (i = 0; i <= deg; i++)
    mpz_fdiv_q_2exp(P[i], P[i], cont);

  return cont;
}

long
bound_roots(mpz_t *t, unsigned long deg)
     /* Johnson's bound : 2*max(abs(-ai/an)^(1/(n-i)),
	the maximum being taken over the i with sgn(a_i) != sgn(an) */
{
  unsigned long i;
  long maxpow, currpow, currpow2, lan, tpos = 1;

  currpow = currpow2 = 0;
  lan = ilog2(t[deg]) - 1; /* puiss de 2 < an */
  maxpow = -lan;

  for (i = 0; i < deg; i++)
    {
      if (mpz_sgn(t[deg]) != mpz_sgn(t[i]))
	{
	  tpos = 0;
	  currpow = ilog2(t[i]);
	  currpow -= lan; /* 2^currpow >= abs(-ai/an) */

	  if (currpow > 0)
	    currpow2 = currpow / (deg - i);
	  else
	    currpow2 = -((-currpow) / (deg - i));

	  if(currpow2 * ((long) (deg - i)) != currpow) currpow2++;
	  /* 2^currpow2 >= abs(-ai/an)^(1/(n-i) */

	  if(currpow2 > maxpow) maxpow = currpow2;
	}
    }

  if (tpos == 1) return -1;

  /* here 2^maxpow > max(abs(-ai/an)^(1/(n-i)), add one to get the bound */
  maxpow++;
  return(maxpow);
}

int
change_pol(mpz_t *P, unsigned long b, unsigned long deg)
     /* From the polynomial P, of degree deg, and the bound b such that
        all positive roots of P are <= 2^b, compute a polynomial Q0
        which has all its real roots in ]0, 1[, namely P(X*2^b). */
{
  long i, j = b;

  for(i = 1; i <= deg; i++, j+=b)
    {
      mpz_mul_2exp(P[i], P[i], j);
    }
  return RemoveContent(P, deg);
}

int
Homoth(mpz_t *P, long k, unsigned long deg)
     /* Computes P(X*2^k) */
{
  long i, j;

#ifdef DEBUG
  fprintf(stderr, "Homothethy of k = %d\n", k);
#endif

  if (k > 0) {
    j = k;
    for(i = 1; i <= deg; i++, j+=k)
      mpz_mul_2exp(P[i], P[i], j);
  }
  else
    {
      j = deg * (-k);
      for(i = 0; i < deg; i++, j += k)
	mpz_mul_2exp(P[i], P[i], j);
    }

  /* Remove possible large power of 2 in content */
  return RemoveContent(P, deg);
}

void
X2XP1(mpz_t *P, unsigned long deg)
     /* Replaces P by the polynomial P(X+1) */
{
  long i, j;

#ifdef DEBUG
  fprintf(stderr, "Translation by 1.\n");
#endif

  for (i = 0; i <= deg-1; i++)
    for (j = deg-1 ; j >= i; j--)
      mpz_add(P[j], P[j], P[j+1]);

  return;
}

/* Number of sign changes in the coefficients of P(1/(X+1))
   (Descartes' rule) */
unsigned long
Descartes(mpz_t *P, unsigned long deg, long sigh, long *flag)
{
  unsigned long nb = 0;
  long i, j, s, t;
  mpz_t *Q;

  node_looked ++;
  /*
     Prune the computation if all the coefficients are of the sign of P[deg]
     In that case any subsequent interval shall have the same property,
     we put *flag at 1 to point this to Uspensky_rec.
  */
  j = deg; t = mpz_sgn(P[j]);
  while (j >= 0 && mpz_sgn(P[j]) == t) { j--; }
  if (j < 0)
    {
#ifdef DEBUG
      fprintf(stderr, "Pruning at i=0, nb = %d\n", nb);
#endif
      *flag = -1;
      return nb;
    }

  Q = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));
  for (i = 0; i <= deg; i++) mpz_init_set(Q[i], P[i]);

  for (j = 0; j <= deg-1; j++)
    mpz_add(Q[j+1], Q[j+1], Q[j]);

  s = mpz_sgn(Q[deg]);

  *flag = s && (s == mpz_sgn(P[0])) && (s == -sigh);

#if DEBUG > 3
  fprintf(stderr, "signe(P[0]) = %d signe(P(1/2)) = %ld signe(P(1)) = %ld, flag = %ld\n", mpz_sgn(P[0]), sigh, s, *flag);
#endif

  for (i = 1; i <= deg-1; i++)
    {
  /*
     Prune the computation if all further coefficients are of the sign of
     Q[deg-i]
  */
      j = deg - i; t = s;
      while (t == 0) { t = mpz_sgn(Q[j]); j--; }
      while (j >= 0 && mpz_sgn(Q[j]) == t) { j--; }
      if (j < 0)
	{
#ifdef DEBUG
	  fprintf(stderr, "Pruning at i=%d, nb = %d\n", i, nb);
#endif
	  for (i = 0; i <= deg; i++) mpz_clear(Q[i]);
	  free(Q);
	  return nb;
	}

      for (j = 0; j <= deg - i - 1; j++)
	mpz_add(Q[j+1], Q[j+1], Q[j]);

      if (s == 0) { s = mpz_sgn(Q[deg-i]); }
      else
	if (s == -mpz_sgn(Q[deg-i]))
	  {
	    if ((nb == 1 && !*flag) || nb == 2)
	      {
		for (i = 0; i <= deg; i++) mpz_clear(Q[i]);
		free(Q);
		return (nb + 1);
	      }

	    nb++; s = -s;
	  }
    }

  if (s == -mpz_sgn(Q[0])) nb++;
  for (i = 0; i <= deg; i++) mpz_clear(Q[i]);

  return nb;
}

/* Returns the sign of P(1/2) */
long
evalhalf(mpz_t *P, unsigned long deg)
{
  long j;
  mpz_t x, y;

  mpz_init_set(x, P[deg]);
  mpz_init(y);

  for (j = deg - 1; j >= 0; j--)
    {
      mpz_mul_2exp(y, P[j], deg - j);
      mpz_add(x, x, y);
    }

  mpz_clear(y);
  return mpz_sgn(x);
}

void
add_root(interval *roots, mpz_t c, int k, unsigned int flag,
	 unsigned int nbroot)
{
  int b = (sign ? b1 : b2);

  mpz_init(roots[nbroot].c);

  if (k <= b)
    {
      if (sign)
	{
	  mpz_neg(roots[nbroot].c, c);
	  if (!flag)
	    mpz_sub_ui(roots[nbroot].c, roots[nbroot].c, 1);
	  mpz_mul_2exp(roots[nbroot].c, roots[nbroot].c, b-k);
	}
      else
	mpz_mul_2exp(roots[nbroot].c, c, b-k);

      roots[nbroot].k = k - b;
      roots[nbroot].isexact = flag;
#ifdef POWER_HACK
      roots[nbroot].sign = 0;
#endif
      return;
    }
  else
    {
      if (sign)
	{
	  mpz_neg(roots[nbroot].c, c);
	  if (!flag)
	    mpz_sub_ui(roots[nbroot].c, roots[nbroot].c, 1);
	}
      else
	mpz_set(roots[nbroot].c, c);

      roots[nbroot].k = k - b;
      roots[nbroot].isexact = flag;
#ifdef POWER_HACK
      roots[nbroot].sign = 0;
#endif
    }
  return;
}

void
affiche_root(FILE *stream, interval z)
{
  mpz_t tmp;

  mpz_init(tmp);
  if (z.isexact != 1) { fprintf(stream, "]"); }

#ifdef POWER_HACK
  if (z.sign != 0 && z.isexact == 0)
    {
      fprintf(stream, "-(");
      if (z.k <= 0)
	{
	  mpz_set_ui(tmp, 1);
	  mpz_mul_2exp(tmp, tmp, -z.k);
	  mpz_add(tmp, z.c, tmp);
	  mpz_out_str(stream, 10, tmp);
	}
      else
	{
	  mpz_add_ui(tmp, z.c, 1);
	  mpz_out_str(stream, 10, tmp); fprintf(stream, "/2^%ld", z.k);
	}
      fprintf(stream, ")^(1/%u)", t_power);
    }
  else
    {
      if (z.sign != 0) fprintf(stream, "-");
      if (t_power > 1) fprintf(stream, "(");
      if (z.k <= 0)
	{ mpz_out_str(stream, 10, z.c); }
      else
	{ mpz_out_str(stream, 10, z.c); fprintf(stream, "/2^%ld", z.k); }
      if (t_power > 1) fprintf(stream, ")^(1/%u)", t_power);
    }
#else

  if (z.k <= 0)
    { mpz_out_str(stream, 10, z.c); }
  else
    { mpz_out_str(stream, 10, z.c); fprintf(stream, "/2^%ld", z.k); }

#endif

  if (z.isexact == 1) { return; }

  fprintf(stream, ", ");


#ifdef POWER_HACK
  if (z.sign != 0)
    {
      fprintf(stream, "-(");
      if (z.k <= 0)
	{ mpz_out_str(stream, 10, z.c); }
      else
	{ mpz_out_str(stream, 10, z.c); fprintf(stream, "/2^%ld", z.k); }
      fprintf(stream, ")^(1/%u)", t_power);
    }
  else
    {
      if (t_power > 1) fprintf(stream, "(");
      if (z.k <= 0)
	{
	  mpz_set_ui(tmp, 1);
	  mpz_mul_2exp(tmp, tmp, -z.k);
	  mpz_add(tmp, z.c, tmp);
	  mpz_out_str(stream, 10, tmp);
	}
      else
	{
	  mpz_add_ui(tmp, z.c, 1);
	  mpz_out_str(stream, 10, tmp); fprintf(stream, "/2^%ld", z.k);
	}
      if (t_power > 1) fprintf(stream, ")^(1/%u)", t_power);
    }
#else
  if (z.k <= 0)
    {
      mpz_set_ui(tmp, 1);
      mpz_mul_2exp(tmp, tmp, -z.k);
      mpz_add(tmp, z.c, tmp);
      mpz_out_str(stream, 10, tmp);
    }
  else
    {
      mpz_add_ui(tmp, z.c, 1);
      mpz_out_str(stream, 10, tmp); fprintf(stream, "/2^%ld", z.k);
    }
#endif

  fprintf(stream, "]");

  mpz_clear(tmp);
}

/*
   Check interval [c/2^k, (c+1)/2^k]. The value of k is returned, this
   is necessary to know from where we come [i.e., what exactly is the
   current polynomial P] when several recursive calls return in a row.
   In practice, this is used to update the polynomial at HERE */

long
Uspensky_rec(mpz_t *P, mpz_t c, unsigned long k,
	     unsigned long *Deg,
	     interval *roots, unsigned int *nbroot)
{
  unsigned long oldk, i, j, nb;
  long shalf, flag;
  mpz_t tmp;

  if (k > max_depth)
    {
      max_depth = k;
      if (k > DEPTH_MAX)
	{
	  fprintf(stderr, "Maximal depth reached. Check that your polynomial is squarefree or increase DEPTH_MAX.\n");
	  exit(-1);
	}
    }

  mpz_init(tmp);

#ifdef DEBUG
  fprintf(stderr, "Checking interval : k = %d, c = ", k);
  mpz_out_str(stderr, 10, c); fprintf(stderr, "\n");
#endif

#if DEBUG > 5
      fprintf(stderr, "Polynomial in this range : ");
      output(P, *Deg);
#endif

      /* Check whether c/2^k is a root */
      if (mpz_cmp_ui(P[0], 0) == 0)
	{
	  i = 1; while(mpz_cmp_ui(P[i], 0) == 0) { i++; }

	  for (j = 0; j < i; j++)
	    {
	      add_root(roots, c, k, 1, *nbroot);
	      (*nbroot) ++;
	    }

	  *Deg -= i; /* Update the polynomial */
	  for (j = 0; j <= *Deg; j++, i++)
	    mpz_set(P[j], P[i]);
	}

      /*
	 Compute the sign of P(1/2) ; thus if Descartes bound is 2,
	 whereas sign(P(0)) = sign(P(1)) = -sign(P(1/2)) we have
	 found two roots.
      */
      shalf = evalhalf(P, *Deg);

      /* Check whether 1/2 is a root */
      while (shalf == 0)
	{
	  /* Print the root. */
	  mpz_set(tmp, c);
	  mpz_mul_2exp(tmp, tmp, 1);
	  mpz_add_ui(tmp, tmp, 1);
	  add_root(roots, tmp, k+1, 1, *nbroot);

	  /* We perform [in place] the division by -1/2 + X */
	  for (i = 0; i < *Deg; i++)
	    {
	      mpz_mul_2exp(P[i], P[i], 1);
	      mpz_add(P[i+1], P[i+1], P[i]);
	      mpz_neg(P[i], P[i]);
	    }
	  (*nbroot)++; (*Deg)--;

	  shalf = evalhalf(P, *Deg);
	}

      /* Compute Descartes' bound */
      nb = Descartes(P, *Deg, shalf, &flag);
      if (flag == TOT_POS)
	{
#ifdef DEBUG
	  fprintf(stderr, "Totally positive, returning.\n");
#endif
	  mpz_clear(tmp);
	  return TOT_POS;
	}

#ifdef DEBUG
      fprintf(stderr, "nb = %ld\n", nb);
#endif

      switch (nb)
	{
	case 0: /* no root */
	  return k;

	case 1: /* exactly one root */
	  add_root(roots, c, k, 0, *nbroot);
	  (*nbroot)++;
	  return k;

	case 2: /* if flag!=0, one root in each half of the current interval */
	  if (flag)
	    {
	      half_done++;
	      mpz_set(tmp, c);

	      mpz_mul_2exp(tmp, tmp, 1);
	      add_root(roots, tmp, k+1, 0, *nbroot);
	      (*nbroot) ++;

	      mpz_add_ui(tmp, tmp, 1);
	      add_root(roots, tmp, k+1, 0, *nbroot);
	      (*nbroot) ++;

	      return k;
	    }

	default: /* recursive call on each half of the interval */
	  mpz_set(tmp, c);

	  mpz_mul_2exp(tmp, tmp, 1);
	  Homoth(P, -1, *Deg);
	  oldk = Uspensky_rec(P, tmp, k+1, Deg, roots, nbroot);
	  if (oldk == TOT_POS)
	    {
#ifdef DEBUG
	      fprintf(stderr, "Totally positive, returning.\n");
#endif
	      mpz_clear(tmp);
	      return TOT_POS;
	    }

	  mpz_add_ui(tmp, tmp, 1);
	  X2XP1(P, *Deg);
	  transl_done++;

	  if (oldk > k + 1)
	    Homoth(P, oldk - (k + 1), *Deg);

	  oldk = Uspensky_rec(P, tmp, k+1, Deg, roots, nbroot);
	  if (oldk == TOT_POS)
	    {
#ifdef DEBUG
	      fprintf(stderr, "Totally positive, returning.\n");
#endif
	      mpz_clear(tmp);
	      return TOT_POS;
	    }

	  mpz_clear(tmp);
	  return oldk;
	}
}

#ifdef POWER_HACK
int gcd
(int a, int b)
{
  int r;

  while (b)
    { r = a%b; a = b; b = r; }
  return a;
}

int
pure_pow(mpz_t *P, unsigned long deg)
{
  int i, t;

  i = deg - 1;
  while (i >= 0 && mpz_sgn(P[i]) == 0) i--;

  t = deg - i;

  for (i--; i >= 0 && t > 1; i--)
    {
      if ((deg - i) % t != 0)
	if (mpz_sgn(P[i])) t = gcd(t, deg - i);
    }

  return t;
}
#endif

interval *
Uspensky(mpz_t *Q, unsigned long deg, unsigned int *nbroot)
{
  interval *roots = (interval *)malloc(deg * sizeof(interval));
  unsigned long deg0 = deg;
  int i;
  mpz_t *P, e;
#ifdef POWER_HACK
  int nb_z;
#endif

  mpz_init_set_ui(e, 0);
  *nbroot = 0;

#ifdef POWER_HACK
  t_power = pure_pow(Q, deg);

  nb_z = deg % t_power;
  while (mpz_sgn(Q[nb_z]) == 0) nb_z += t_power;


  for (j = 0; j < nb_z; j++)
    {
      add_root(roots, e, 0, 1, *nbroot);
      (*nbroot) ++;
    }

  deg = (deg0 - nb_z)/t_power; /* Update the polynomial */
  P = (mpz_t *) malloc ((deg + 1) * sizeof(mpz_t));

  for (j = 0; j <= deg; j++)
    mpz_init_set(P[j], Q[nb_z + t_power*j]);

#ifdef DEBUG
  fprintf(stderr, "POWER_HACK: power = %d\n", t_power);
  fprintf(stderr, "POWER_HACK: 0 root of multiplicity %d\n", nb_z);
  fprintf(stderr, "POWER_HACK: polynomial is now "); output(P, deg);
#endif
#else /* ifndef POWER_HACK */
  P = (mpz_t *) malloc ((deg + 1) * sizeof(mpz_t));

  for (i = deg; i >= 0; i--)
    mpz_init_set(P[i], Q[i]);
#endif

  /* First work on the positive roots. */
  if (!b_pos_flag)
    b2 = bound_roots(P, deg);

  if (bound_pos_flag)
    {
      if (b2 >= 0)
	{
	  fprintf(stderr, "Using bound 2^%d on positive roots", b2);
#ifdef POWER_HACK
	  if (t_power > 1) fprintf(stderr, " (w. power_hack)");
#endif
	  fprintf(stderr, ".\n");
	}
      else
	{
          fprintf(stderr, "Using bound 0 on positive roots ");
#ifdef POWER_HACK
	  if (t_power > 1) fprintf(stderr, " (w. power_hack)");
#endif
	  fprintf(stderr, ".\n");
	}
    }

  if (b2 < 0) goto NEGATIVE;

  change_pol(P, b2, deg);

  sign = 0;
  Uspensky_rec(P, e, 0, &deg, roots, nbroot);

#ifdef POWER_HACK /* if t is even, add opposites of nonzero roots */
  if ((t_power & 1) == 0)
    {
      for (i = 0; i < *nbroot - nb_z; i++)
	{
	  roots[*nbroot + i].k = roots[i + nb_z].k;
	  mpz_init_set(roots[*nbroot + i].c, roots[i + nb_z].c);
	  roots[*nbroot + i].isexact = roots[i + nb_z].isexact;
	  roots[*nbroot + i].sign = 1;
	}
      (*nbroot) = 2 * (*nbroot) - nb_z;
    }

  /* Change P into P(-X) to look for negative roots */
 NEGATIVE:
  deg = (deg0 - nb_z) /t_power;
  for (j = 0; j <= deg; j++)
    {
      if (j % 2 == 1)
	mpz_neg(P[j], Q[nb_z + t_power*j]);
      else mpz_set(P[j], Q[nb_z + t_power*j]);
    }
#else
  /* Change P into P(-X) to look for negative roots */
 NEGATIVE:
  deg = deg0;
  for (i = deg; i >= 0; i--)
    {
      if (i % 2 == 1)
	mpz_neg(P[i], Q[i]);
      else mpz_set(P[i], Q[i]);
    }
#endif

  if (!b_neg_flag)
    b1 = bound_roots(P, deg);

  if (bound_neg_flag)
    {
      if (b1 >= 0)
	{
	  fprintf(stderr, "Using bound 2^%d on negative roots", b1);
#ifdef POWER_HACK
	  if (t_power > 1) fprintf(stderr, " (w. power_hack)");
#endif
	  fprintf(stderr, ".\n");
	  change_pol(P, b1, deg);
	  mpz_set_ui(e, 0);
	  sign = 1;
	  Uspensky_rec(P, e, 0, &deg, roots, nbroot);
	}
      else
	{
	  fprintf(stderr, "Using bound 0 on negative roots");
#ifdef POWER_HACK
	  if (t_power > 1) fprintf(stderr, " (w. power_hack)");
#endif
	  fprintf(stderr, "\n");
	}
    }

  /* Free memory. */
#ifndef POWER_HACK
  for (i = deg0; i >= 0; i--)
#else
  for (i = (deg0 - nb_z) / t_power; i >= 0; i--)
#endif
    mpz_clear(P[i]);

  free(P);

  return roots;
}

#if 0

int
main(int argc, char **argv)
{
  mpz_t *P;
  interval *roots;
  int i;
  unsigned long deg;
  unsigned int nbroot, nb_flag = 0, rev_flag = 0;
  unsigned int roots_flag = 1, stats_flag = 0;
  FILE *file_in = stdin;
  char *progname = argv[0];

  argc--; argv++;
  while (argc) {
    if (!strcmp(*argv, "-n"))
      { argv++; argc--; nb_flag = 1; }
    else if (!strcmp(*argv, "-d"))
      { argv++; argc--; rev_flag = 1; }
    else if (!strcmp(*argv, "-a"))
      { argv++; argc--; rev_flag = 0; }
    else if (!strcmp(*argv, "-x"))
      { argv++; argc--; roots_flag = 0; }
    else if (!strcmp(*argv, "-b"))
      {
	argv++; argc--;
	bound_neg_flag = bound_pos_flag = 1;
      }
    else if (!strcmp(*argv, "-b1"))
      {
	argv++; argc--;
	bound_neg_flag = 1;
	if (argc && argv[0][0] != '-')
	  {
	    b_neg_flag = 1;
	    b1 = atoi(argv[0]);
	    argv++; argc--;
	  }
      }
    else if (!strcmp(*argv, "-b2"))
      {
	argv++; argc--;
	bound_pos_flag = 1;
	if (argc && argv[0][0] != '-')
	  {
	    b_pos_flag = 1;
	    b2 = atoi(argv[0]);
	    argv++; argc--;
	  }
      }
    else if (!strcmp(*argv, "-s"))
      {
	argv++; argc--; stats_flag = 1;
	transl_done = node_looked = half_done = max_depth = 0;
      }
    else if (!strcmp(*argv, "-f"))
      {
	argv++; argc--;
	if ((file_in = fopen(argv[0], "r")) == NULL)
	  {
	    fprintf(stderr, "Unable to open file %s.\n", argv[0]);
	    exit(-1);
	  }
	argv++; argc--;
      }
    else {
      if (strcmp(*argv, "-h"))
	fprintf(stderr,"Unknown option : %s\n",*argv);

      fprintf(stderr,"syntax is : %s [-hndabsx] [-f file] [-b1 int] [-b2 int]\n", progname);
      fprintf(stderr,"-h: print this message\n");
      fprintf(stderr,"-d: read coefficients in descending order\n");
      fprintf(stderr,"-a: read coefficients in ascending order [default]\n");
      fprintf(stderr,"-f file: read coefficients from file [default stdin]\n");
      fprintf(stderr,"-n: print the number of roots [default : off]\n");
      fprintf(stderr,"-s: print various statistics [default : off]\n");
      fprintf(stderr,"-x: do not print the roots [default : off]\n");
      fprintf(stderr,"-b: print the bounds on positive and negative roots.\n");
      fprintf(stderr,"-b1 [int]: print the bound on negative roots, or use the given nonnegative integer as a bound.\n");
      fprintf(stderr,"-b2 [int]: print the bound on positive roots, or use the given nonnegative integer as a bound.\n");
      exit(0);
    }
  }

  fscanf(file_in, "%lu\n", &deg);
  P = (mpz_t *)malloc ((deg + 1) * sizeof(mpz_t));

  if (rev_flag)
    for (i = deg; i >= 0; i--)
      {
	mpz_init(P[i]);
	mpz_inp_str(P[i], file_in, 10);
      }
  else
    for (i = 0; i <= deg; i++)
      {
	mpz_init(P[i]);
	mpz_inp_str(P[i], file_in, 10);
      }

  roots = Uspensky(P, deg, &nbroot);

  if (nb_flag)
    fprintf(stdout, "%u roots\n", nbroot);

  if (roots_flag)
    {
      for (i = 0; i < nbroot; i++)
	{
	  affiche_root(stdout, roots[i]);
	  if (i < nbroot - 1) fprintf(stdout, ", ");
	}
      fprintf(stdout, "\n");
    }

  if (stats_flag)
    {
      fprintf(stderr, "Nodes visited = %lu\n", node_looked);
      fprintf(stderr, "Max depth = %lu\n", max_depth);
      fprintf(stderr, "Translations performed = %lu\n", transl_done);
      fprintf(stderr, "Nodes with 1/2 hack = %lu\n", half_done);
    }
  return(0);
}

#endif
