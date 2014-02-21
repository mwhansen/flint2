/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"


/* Computes

   xi    = x ^ (p^(i*s)),
   zeta  = lambda ^ (p^s + p^(2*s) + ... + p^(i*s))
   delta =
   
   using Algorithm 1 (XiZetaDelta).
*/

void
__TEMPLATE(T, root_ds_xizetadelta)(TEMPLATE(T, t) xi, TEMPLATE(T, t) zeta,
                                      TEMPLATE(T, t) delta,
                                      const TEMPLATE(T, t) op,
                                      ulong i,
                                      const TEMPLATE(T, t) xi1,
                                      const TEMPLATE(T, t) zeta1,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    int j;
    TEMPLATE(T, t) t1, t2;

    TEMPLATE(T, set)(xi, xi1, ctx);
    TEMPLATE(T, set)(zeta, zeta1, ctx);
    TEMPLATE(T, set)(delta, zeta1, ctx);

    TEMPLATE(T, init)(t1, ctx);
    TEMPLATE(T, init)(t2, ctx);
    for (j = ((int)FLINT_BIT_COUNT(i) - 2); j >= 0; j--)
    {
        /* delta += zeta * delta(xi) */
        TEMPLATE(T, compose)(t1, delta, xi, ctx);
        TEMPLATE(T, mul)(t2, zeta, t1, ctx);
        TEMPLATE(T, add)(delta, delta, t2, ctx);

        /* zeta = zeta * zeta(xi) */
        TEMPLATE(T, compose)(t1, zeta, xi, ctx);
        TEMPLATE(T, mul)(t2, zeta, t1, ctx);
        TEMPLATE(T, swap)(t2, zeta, ctx);

        /* xi = xi(xi) */
        TEMPLATE(T, compose)(t1, xi, xi, ctx);
        TEMPLATE(T, swap)(xi, t1, ctx);

        if (i & UWORD(1) << j)
        {
            TEMPLATE(T, compose)(t1, xi, xi1, ctx);
            TEMPLATE(T, swap)(t1, xi, ctx);

            TEMPLATE(T, compose)(t1, zeta, xi1, ctx);
            TEMPLATE(T, mul)(t2, zeta1, t1, ctx);
            TEMPLATE(T, swap)(zeta, t2, ctx);

            TEMPLATE(T, add)(delta, delta, zeta, ctx);
        }
        
    }
    
    TEMPLATE(T, clear)(t1, ctx);
    TEMPLATE(T, clear)(t2, ctx);
}


/* Computes

   xi   = x ^ (p^(i*s)),
   zeta = lambda ^ (p^s + p^(2*s) + ... + p^(i*s))

   using Algorithm 1 (XiZetaDelta) (without the delta).
*/
void
__TEMPLATE(T, root_ds_xizeta)(TEMPLATE(T, t) xi, TEMPLATE(T, t) zeta,
                              ulong i,
                              const TEMPLATE(T, t) xi1,
                              const TEMPLATE(T, t) zeta1,
                              const TEMPLATE(T, ctx_t) ctx)
{
    int j;
    TEMPLATE(T, t) t1, t2;

    TEMPLATE(T, set)(xi, xi1, ctx);
    TEMPLATE(T, set)(zeta, zeta1, ctx);

    TEMPLATE(T, init)(t1, ctx);
    TEMPLATE(T, init)(t2, ctx);
    for (j = ((int)FLINT_BIT_COUNT(i) - 2); j >= 0; j--)
    {
        /* zeta = zeta * zeta(xi) */
        TEMPLATE(T, compose)(t1, zeta, xi, ctx);
        TEMPLATE(T, mul)(t2, zeta, t1, ctx);
        TEMPLATE(T, swap)(zeta, t2, ctx);

        /* xi = xi(xi) */
        TEMPLATE(T, compose)(t1, xi, xi, ctx);
        TEMPLATE(T, swap)(xi, t1, ctx);

        if (i & UWORD(1) << j)
        {
            /* xi = xi(xi1) */
            TEMPLATE(T, compose)(t1, xi, xi1, ctx);
            TEMPLATE(T, swap)(xi, t1, ctx);

            /* zeta = zeta1 * zeta(xi1) */
            TEMPLATE(T, compose)(t1, zeta, xi1, ctx);
            TEMPLATE(T, mul)(t2, zeta1, t1, ctx);
            TEMPLATE(T, swap)(zeta, t2, ctx);
        }
    }
    
    TEMPLATE(T, clear)(t1, ctx);
    TEMPLATE(T, clear)(t2, ctx);
}

/* Computes xp    = x ^ p
            xi1   = x ^ (p^s)
            zeta1 = op ^ (p ^ s)
 */
void
__TEMPLATE(T, root_ds_xi1zeta1)(TEMPLATE(T, t) xi1, TEMPLATE(T, t) zeta1,
                                TEMPLATE(T, t) xp,
                                const TEMPLATE(T, t) op,
                                const fmpz_t s,
                                const TEMPLATE(T, ctx_t) ctx)
{
    int j;
    
    /* xp = x ^ p */
    TEMPLATE(T, gen)(xi1, ctx);
    TEMPLATE(T, pow)(xp, xi1, TEMPLATE(T, ctx_prime)(ctx), ctx);

    /* xi1 = x ^ p */
    TEMPLATE(T, set)(xi1, xp, ctx);

    /* xi1 = x ^ (p ^ s) */
    for (j = fmpz_sizeinbase(s, 2) - 2; j >= 0; j--)
    {
        TEMPLATE(T, compose)(xi1, xi1, xi1, ctx);
        if (fmpz_tstbit(s, j))
            TEMPLATE(T, compose)(xi1, xi1, xp, ctx);
    }

    TEMPLATE(T, compose)(zeta1, op, xi1, ctx);
}


void
__TEMPLATE(T, root_ds_alpha)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                             ulong i, fmpz_t s, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) xi1, zeta1, xp, xi, zeta, delta;

    TEMPLATE(T, init)(xi1, ctx);
    TEMPLATE(T, init)(zeta1, ctx);
    TEMPLATE(T, init)(xp, ctx);
    __TEMPLATE(T, root_ds_xi1zeta1)(xi1, zeta1, xp, op, s, ctx);
    
    TEMPLATE(T, init)(xi, ctx);
    TEMPLATE(T, init)(zeta, ctx);
    TEMPLATE(T, init)(delta, ctx);

    __TEMPLATE(T, root_ds_xizetadelta)(xi, zeta, delta, op, i,
                                       xi1, zeta1, ctx);

    TEMPLATE(T, mul)(rop, op, delta, ctx);

    TEMPLATE(T, clear)(xi, ctx);
    TEMPLATE(T, clear)(zeta, ctx);
    TEMPLATE(T, clear)(xi1, ctx);
    TEMPLATE(T, clear)(zeta1, ctx);
    TEMPLATE(T, clear)(delta, ctx);
    TEMPLATE(T, clear)(xp, ctx);
}

void
__TEMPLATE(T, root_ds_s_ell)(fmpz_t s, fmpz_t ell, fmpz_t e, slong t,
                             const TEMPLATE(T, ctx_t) ctx)
{
    fmpz_t ft, p, n;
    fmpz_init(ft);
    fmpz_init(p);
    fmpz_init_set_ui(n, TEMPLATE(T, ctx_degree)(ctx));
    
    fmpz_set_si(ft, t);
    fmpz_mod(p, TEMPLATE(T, ctx_prime)(ctx), ft);
    fmpz_gcd(ell, ft, p);
    fmpz_divexact(s, n, ell);

    fmpz_pow_ui(e, TEMPLATE(T, ctx_prime)(ctx), fmpz_get_ui(s));
    fmpz_sub_ui(e, e, 1);
    fmpz_divexact_si(e, e, t);

    fmpz_clear(ft);
    fmpz_clear(p);
    fmpz_clear(n);
}


/* Computes op ^ ((q - 1) / t) */
void
__TEMPLATE(T, root_ds_qm1ot)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                             const fmpz_t s, const fmpz_t ell, const fmpz_t e,
                             slong t, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) zeta, xp, xi, xi1, zeta1;

    /* Handle the case where ell == 1 <=> p^s == q */
    /* This shouldn't happen if we call root_edf when s == d */
    if (fmpz_cmp_ui(ell, 1) == 0)
    {
        TEMPLATE(T, pow)(rop, op, e, ctx);
        return;
    }
    
    /* TODO: Profile to find cutoff vs. direct powering */
    
    /* e = (p^s - 1) / t
       zeta = op ^ ((p^s - 1) / t) */
    TEMPLATE(T, init)(zeta, ctx);
    TEMPLATE(T, pow)(zeta, op, e, ctx);

    /* Use Algorithm 1 to compute
       op ^ ((q - 1) / t) = zeta^(1 + p^s + ... + p^(s*(l - 1))) */
    TEMPLATE(T, init)(xp, ctx);
    TEMPLATE(T, init)(xi1, ctx);
    TEMPLATE(T, init)(zeta1, ctx);

    /* Compute xi1   = x ^ (p^s)
               zeta1 = zeta ^ (p^s) */ 
    __TEMPLATE(T, root_ds_xi1zeta1)(xi1, zeta1, xp, zeta, s, ctx);

    /* Use a variant of Algorithm 1 (XiZetaDelta) to compute
       rop = zeta ^ (p^s + p^(2*s) + ... + p^(i*s)) */
    TEMPLATE(T, init)(xi, ctx);
    __TEMPLATE(T, root_ds_xizeta)(xi, rop, fmpz_get_ui(ell) - 1, xi1, zeta1,
                                  ctx);

    /* rop = zeta^(1 + p^s + p^(2*s) + ... + p^(i*s))
           = op ^ (q - 1) / t */
    TEMPLATE(T, mul)(rop, rop, zeta, ctx);

    TEMPLATE(T, clear)(zeta, ctx);
}

/* Assumes that op is a t-th power */
void
_TEMPLATE(T, root_ds)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                      slong t, const TEMPLATE(T, ctx_t) ctx)

{
    int res;
    flint_rand_t state;
    fmpz_t e, s, ell, ft;
    TEMPLATE(T, ctx_t) ctxp;
    TEMPLATE(T, t) ap, apbt, b, c, t1, lambda, beta, z;

    fmpz_init(s);
    fmpz_init(e);
    fmpz_init(ell);
    fmpz_init(ft);
    fmpz_set_si(ft, t);
    __TEMPLATE(T, root_ds_s_ell)(s, ell, e, t, ctx);
    
    flint_randinit(state);
    TEMPLATE(T, init)(ap, ctx);
    TEMPLATE(T, init)(b, ctx);
    TEMPLATE(T, init)(c, ctx);
    TEMPLATE(T, init)(apbt, ctx);
    TEMPLATE(T, init)(t1, ctx);
    TEMPLATE(T, init)(lambda, ctx);
    while (TEMPLATE(T, is_zero)(b, ctx))
    {
        /* Random c */
        TEMPLATE(T, randtest)(c, state, ctx);

        /* ap = op * c ^ t */
        TEMPLATE(T, pow)(t1, c, ft, ctx);
        TEMPLATE(T, mul)(ap, op, t1, ctx);

        /* lambda = ap ^ ((p^s - 1) / t) */
        TEMPLATE(T, pow)(lambda, ap, e, ctx);

        /* b = 1 + lambda + alpha(lambda, ell - 2) */
        TEMPLATE(T, one)(b, ctx);
        TEMPLATE(T, add)(b, b, lambda, ctx);
        /*TODO: is ell always >= 2*/
        __TEMPLATE(T, root_ds_alpha)(t1, lambda, fmpz_get_ui(ell) - 2, s, ctx);
        TEMPLATE(T, add)(b, b, t1, ctx);
    }

    TEMPLATE(T, pow_ui)(t1, b, t, ctx);
    TEMPLATE(T, mul)(apbt, ap, t1, ctx);

    /* Check to see if we can just work over the prime field */
    if (TEMPLATE(T, is_in_prime_field)(apbt, ctx))
    {
        flint_printf("Prime field case\n");
        _TEMPLATE(T, root_prime)(rop, apbt, t, ctx);
    }
    else
    {
        TEMPLATE(T, ctx_init_minpoly)(ctxp, apbt, "z", ctx);
        TEMPLATE(T, init)(z, ctxp);
        TEMPLATE(T, init)(beta, ctxp);
        TEMPLATE(T, gen)(z, ctxp);

        /* Compute the t^th root in the smaller field */
        res = TEMPLATE(T, root_edf)(beta, z, t, ctxp);
        flint_printf("Subfield case: %d\n", res);
        TEMPLATE(T, ctx_print)(ctx);
        flint_printf("d = %wd\ns = ", TEMPLATE(T, ctx_degree)(ctx)); fmpz_print(s);
        flint_printf("\n");

        if (res == 0)
        {
            TEMPLATE(T, ctx_print)(ctxp);
            TEMPLATE(T, print_pretty)(z, ctx);
            flint_printf("\n");
        }

        /* rop = Embed(beta, ctx) */
        TEMPLATE(T, compose)(rop, beta, apbt, ctx);

        TEMPLATE(T, clear)(z, ctxp);
        TEMPLATE(T, clear)(beta, ctxp);    
        TEMPLATE(T, ctx_clear)(ctxp);
    }

    /* rop = rop / (b * c) */
    /* We don't need lambda anymore so we can use it as a temporary */
    TEMPLATE(T, mul)(lambda, b, c, ctx);
    TEMPLATE(T, inv)(t1, lambda, ctx);
    TEMPLATE(T, mul)(lambda, t1, rop, ctx);
    TEMPLATE(T, swap)(rop, lambda, ctx);


    TEMPLATE(T, clear)(ap, ctx);
    TEMPLATE(T, clear)(b, ctx);
    TEMPLATE(T, clear)(c, ctx);
    TEMPLATE(T, clear)(apbt, ctx);
    TEMPLATE(T, clear)(t1, ctx);
    TEMPLATE(T, clear)(lambda, ctx);
    
    fmpz_clear(s);
    fmpz_clear(ell);
    fmpz_clear(ft);
    flint_randclear(state);
}                       

int
TEMPLATE(T, root_ds)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                        slong t, const TEMPLATE(T, ctx_t) ctx)
{
    int res = 0;
    fmpz_t qm1ot, s, ell, e;
    TEMPLATE(T, t) at;

    /* Test if op is indeed a t-th power */
    /* TODO: Make this much more efficient when t | p - 1*/
    fmpz_init(qm1ot);
    TEMPLATE(T, init)(at, ctx);
    
    TEMPLATE(T, ctx_order)(qm1ot, ctx);
    fmpz_sub_ui(qm1ot, qm1ot, 1);
    if (!fmpz_divisible_si(qm1ot, t)) /* Check if t | q - 1 */
    {
        res = 0;
        goto cleanup;
    }
    else
    {
        fmpz_divexact_si(qm1ot, qm1ot, t);
    }

    /* Compute s, ell, e */
    fmpz_init(s);
    fmpz_init(ell);
    fmpz_init(e); /* e = (p^s - 1) / t */
    __TEMPLATE(T, root_ds_s_ell)(s, ell, e, t, ctx);

    if (fmpz_cmp_ui(s, TEMPLATE(T, ctx_degree)(ctx)) == 0)
    {
        res = TEMPLATE(T, root_edf)(rop, op, t, ctx);
        goto cleanup;
    }

    /* Compute at = op ^ (q - 1) / t */
    __TEMPLATE(T, root_ds_qm1ot)(at, op, s, ell, e, t, ctx);

    if (!TEMPLATE(T, is_one)(at, ctx))
    {
        res = 0;
        goto cleanup;
    }
    else
    {
        res = 1;
    }

    /* Make the real call now that we know that t | q - 1 and op is
       indeed a t-th power */
    _TEMPLATE(T, root_ds)(rop, op, t, ctx);

cleanup:
    fmpz_clear(qm1ot);
    TEMPLATE(T, clear)(at, ctx);
    
    return res;

}

#endif
