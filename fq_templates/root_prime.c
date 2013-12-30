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


void
__TEMPLATE(T, root_prime_xizetadelta)(TEMPLATE(T, t) xi, TEMPLATE(T, t) zeta,
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

void
__TEMPLATE(T, root_prime_alpha)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                                ulong i, fmpz_t s, const TEMPLATE(T, ctx_t) ctx)
{
    int j;
    TEMPLATE(T, t) xi1, zeta1, xp, xi, zeta, delta;

    TEMPLATE(T, init)(xi1, ctx);
    TEMPLATE(T, init)(zeta1, ctx);
    TEMPLATE(T, init)(xp, ctx);

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
            TEMPLATE(T, compose(xi1, xi1, xp, ctx));
    }

    TEMPLATE(T, compose)(zeta1, op, xi1, ctx);

    TEMPLATE(T, init)(xi, ctx);
    TEMPLATE(T, init)(zeta, ctx);
    TEMPLATE(T, init)(delta, ctx);

    __TEMPLATE(T, root_prime_xizetadelta)(xi, zeta, delta, op, i,
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
TEMPLATE(T, root_prime)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                        ulong t, const TEMPLATE(T, ctx_t) ctx)

{
    flint_rand_t state;
    fmpz_t s, ell, ft;
    TEMPLATE(T, ctx_t) ctxp;
    TEMPLATE(T, t) zeta, ap, apbt, b, c, t1, lambda, beta, z;

    fmpz_init(s);
    fmpz_init(ell);
    fmpz_init_set_ui(ft, t);
    fmpz_gcd(ell, ft, TEMPLATE(T, ctx_prime)(ctx));
    fmpz_divexact(s, ft, ell);

    /* Test if op is indeed a t-th root */
    
    /* zeta = op ^ (p^s - 1) */

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

        /* b = 1 + lambda + alpha(lambda, ell - 2) */
        TEMPLATE(T, one)(b, ctx);
        TEMPLATE(T, add)(b, b, lambda, ctx);
        __TEMPLATE(T, root_prime_alpha)(t1, lambda, fmpz_get_ui(ell) - 2, s, ctx);
        TEMPLATE(T, add)(b, b, t1, ctx);
    }

    TEMPLATE(T, pow_ui)(t1, b, t, ctx);
    TEMPLATE(T, mul)(apbt, ap, t1, ctx); /* We don't need lambda anymore */
    TEMPLATE(T, ctx_init_minpoly)(ctxp, apbt, "z", ctx);

    TEMPLATE(T, init)(z, ctxp);
    TEMPLATE(T, init)(beta, ctxp);
    TEMPLATE(T, one)(z, ctxp);
    TEMPLATE(T, root_edf)(beta, z, t, ctxp);

    /* rop = Embed(beta, ctx) * b^1 * c^-1 */
    TEMPLATE(T, compose)(rop, beta, apbt, ctx);
    TEMPLATE(T, mul)(lambda, b, c, ctx);
    TEMPLATE(T, inv)(t1, lambda, ctx);
    TEMPLATE(T, mul)(lambda, t1, rop, ctx);
    TEMPLATE(T, swap)(lambda, rop, ctx);
    
    
    TEMPLATE(T, ctx_clear)(ctxp);
    fmpz_clear(s);
    fmpz_clear(ell);
    fmpz_clear(ft);
    flint_randclear(state);
}                       

#endif
