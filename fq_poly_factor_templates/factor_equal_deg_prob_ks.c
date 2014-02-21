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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"
#include "ulong_extras.h"

/* Algorithm T from "Fast Polynomial Factorization Over High Algebraic
   Extensions of Finite Fields" by Kaltofen and Shoup. */
void
__TEMPLATE(T, poly_factor_equal_deg_prob_ks_t)(TEMPLATE(T, poly_t) beta,
                                               TEMPLATE(T, poly_t) xi,
                                               TEMPLATE(T, t) zeta,
                                               slong i,
                                               const TEMPLATE(T, poly_t) alpha,
                                               const TEMPLATE(T, poly_t) f,
                                               const TEMPLATE(T, poly_t) finv,
                                               const TEMPLATE(T, ctx_t) ctx)
{
    int j;
    slong ell;
    TEMPLATE(T, poly_t) beta_1, xi_1, temp1, temp2;
    TEMPLATE(T, t) zeta_1;
    
    
    /* beta_1 = alpha ^ p */
    TEMPLATE(T, poly_init)(beta_1, ctx);
    TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(beta_1, alpha,
                                                 TEMPLATE(T, ctx_prime)(ctx),
                                                 0, f, finv, ctx);
    TEMPLATE(T, poly_set)(beta, beta_1, ctx);

    /* xi_1 = x ^ p mod f */
    TEMPLATE(T, poly_init)(xi_1, ctx);
    TEMPLATE(T, poly_gen)(xi, ctx);
    TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(xi_1, xi,
                                                 TEMPLATE(T, ctx_prime)(ctx),
                                                 0, f, finv, ctx);
    TEMPLATE(T, poly_set)(xi, xi_1, ctx);
    
    /* zeta_1 = a ^ p */
    TEMPLATE(T, init)(zeta_1, ctx);
    TEMPLATE(T, gen)(zeta, ctx);
    TEMPLATE(T, pow)(zeta_1, zeta, TEMPLATE(T, ctx_prime)(ctx), ctx);
    TEMPLATE(T, set)(zeta, zeta_1, ctx);
    
    if (i == 1)
        goto cleanup;

    TEMPLATE(T, poly_init)(temp1, ctx);
    TEMPLATE(T, poly_init)(temp2, ctx);
    for (j = ((int)FLINT_BIT_COUNT(i) - 2); j >= 0; j--)
    {
        /* Step T2 (temp1 = beta ^ (p^j)) */
        TEMPLATE(T, poly_fit_length)(temp1, beta->length, ctx);
        for (ell = 0; ell < beta->length; ell++)
        {
            TEMPLATE(T, compose)(temp1->coeffs + ell,
                                 beta->coeffs + ell,
                                 zeta, ctx);
        }
        _TEMPLATE(T, poly_set_length)(temp1, beta->length, ctx);
        _TEMPLATE(T, poly_normalise)(temp1, ctx);
        TEMPLATE(T, poly_compose_mod_preinv)(temp2, temp1, xi, f, finv, ctx);
        TEMPLATE(T, poly_add)(beta, beta, temp2, ctx);

        /* Step T2 (xi = xi ^ (p^j)) */
        TEMPLATE(T, poly_fit_length)(temp2, xi->length, ctx);
        for (ell = 0; ell < xi->length; ell++)
        {
            TEMPLATE(T, compose)(temp2->coeffs + ell,
                                 xi->coeffs + ell,
                                 zeta, ctx);
        }
        _TEMPLATE(T, poly_set_length)(temp2, xi->length, ctx);
        _TEMPLATE(T, poly_normalise)(temp2, ctx);
        TEMPLATE(T, poly_compose_mod_preinv)(xi, temp2, xi, f, finv, ctx);

        /* Step T2 (zeta = zeta ^ (p^j)) */
        TEMPLATE(T, compose)(zeta, zeta, zeta, ctx);

        /* Step T3 */
        if (i & (WORD(1) << j))
        {
            /* Compute beta_{2*j + 1} */
            TEMPLATE(T, poly_fit_length)(temp1, beta->length, ctx);
            for (ell = 0; ell < beta->length; ell++)
            {
                TEMPLATE(T, compose)(temp1->coeffs + ell,
                                     beta->coeffs + ell,
                                     zeta_1, ctx);
            }
            _TEMPLATE(T, poly_set_length)(temp1, beta->length, ctx);
            _TEMPLATE(T, poly_normalise)(temp1, ctx);
            TEMPLATE(T, poly_compose_mod_preinv)(temp2, temp1, xi_1, f, finv, ctx);
            TEMPLATE(T, poly_add)(beta, beta_1, temp2, ctx);

            /* Compute xi = xi ^ p */
            TEMPLATE(T, poly_fit_length)(temp2, xi->length, ctx);
            for (ell = 0; ell < xi->length; ell++)
            {
                TEMPLATE(T, compose)(temp2->coeffs + ell,
                                     xi->coeffs + ell,
                                     zeta_1, ctx);
            }
            _TEMPLATE(T, poly_set_length)(temp2, xi->length, ctx);
            _TEMPLATE(T, poly_normalise)(temp2, ctx);
            TEMPLATE(T, poly_compose_mod_preinv)(xi, temp2, xi_1, f, finv, ctx);

            /* Compute zeta = zeta ^ p */
            TEMPLATE(T, compose)(zeta, zeta, zeta_1, ctx);
        }
    }

    TEMPLATE(T, poly_clear)(temp1, ctx);
    TEMPLATE(T, poly_clear)(temp2, ctx);

cleanup:
    TEMPLATE(T, poly_clear)(xi_1, ctx);
    TEMPLATE(T, poly_clear)(beta_1, ctx);
    TEMPLATE(T, clear)(zeta_1, ctx);
}

/* Algorithm E from "Fast Polynomial Factorization Over High Algebraic
Extensions of Finite Fields" by Kaltofen and Shoup. */

int
TEMPLATE(T, poly_factor_equal_deg_prob_ks) (TEMPLATE(T, poly_t) f,
                                            TEMPLATE(T, poly_t) g,
                                            flint_rand_t state,
                                            const TEMPLATE(T, poly_t) pol,
                                            slong d,
                                            const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_struct) *h2;
    TEMPLATE(T, poly_t) a, gamma, polinv, beta, xi;
    TEMPLATE(T, t) zeta;
    fmpz_t exp;
    int res = 0;

    if (pol->length <= 1)
    {
        TEMPLATE_PRINTF("Exception (%s_poly_factor_equal_deg_prob_ks): \n", T);
        flint_printf("Input polynomial is linear.\n");
        abort();
    }

    /* Compute random $a = \alpha \in \mathbb{A}_{pol}$ */
    TEMPLATE(T, poly_init) (a, ctx);
    do
    {
        TEMPLATE(T, poly_randtest) (a, state, pol->length - 1, ctx);
    } while (a->length <= 1);

    TEMPLATE(T, poly_init) (polinv, ctx);
    TEMPLATE(T, poly_reverse) (polinv, pol, pol->length, ctx);
    TEMPLATE(T, poly_inv_series_newton) (polinv, polinv, polinv->length, ctx);

    
    /* Step E1: Compute trace-like map
       a = a + a^p + a^(p^2) + ... + a^(p^(k*d - 1)) */
    TEMPLATE(T, poly_init) (beta, ctx);
    TEMPLATE(T, poly_init) (xi, ctx);
    TEMPLATE(T, init) (zeta, ctx);
    __TEMPLATE(T, poly_factor_equal_deg_prob_ks_t)(beta, xi, zeta,
                                                   TEMPLATE(T, ctx_degree)(ctx)*d - 1,
                                                   a, pol, polinv, ctx);
    TEMPLATE(T, poly_add)(beta, beta, a, ctx);

    /* Step E2: */
    /* If p > 2, compute beta^((p - 1)/2) */
    TEMPLATE(T, poly_init) (gamma, ctx);
    if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime) (ctx), 2) > 0)
    {
        /* compute a^{(p-1)/2} rem pol */
        fmpz_init(exp);
        fmpz_sub_ui(exp, TEMPLATE(T, ctx_prime)(ctx), 1);
        fmpz_fdiv_q_2exp(exp, exp, 1);

        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (gamma, beta, exp, 0, pol,
                                                      polinv, ctx);
        fmpz_clear(exp);
    }
    else
    {
        TEMPLATE(T, poly_set)(gamma, beta, ctx);
    }

    TEMPLATE(T, poly_gcd)(f, gamma, pol, ctx);
    if (!TEMPLATE(T, poly_is_one)(f, ctx) &&
        !TEMPLATE(T, poly_equal)(f, pol, ctx))
    {
        h2 = g;
        res += 1;
    }
    else
    {
        h2 = f;
    }

    TEMPLATE(T, poly_one)(xi, ctx);
    TEMPLATE(T, poly_add)(gamma, gamma, xi, ctx);
    TEMPLATE(T, poly_gcd)(h2, gamma, pol, ctx);
    if (!TEMPLATE(T, poly_is_one)(h2, ctx) && 
        !TEMPLATE(T, poly_equal)(h2, pol, ctx))
    {
        res += 1;
    }

    TEMPLATE(T, poly_clear) (a, ctx);
    TEMPLATE(T, poly_clear) (beta, ctx);
    TEMPLATE(T, poly_clear) (xi, ctx);
    TEMPLATE(T, poly_clear) (gamma, ctx);
    TEMPLATE(T, clear) (zeta, ctx);
    TEMPLATE(T, poly_clear) (polinv, ctx);

    return res;
}


#endif
