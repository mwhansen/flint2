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

int
TEMPLATE(T, poly_factor_equal_deg_prob_vzgs) (TEMPLATE(T, poly_t) f,
                                              TEMPLATE(T, poly_t) g,
                                              flint_rand_t state,
                                              const TEMPLATE(T, poly_t) pol,
                                              slong d,
                                              const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_struct) *h2;
    TEMPLATE(T, poly_t) a, b, c, xq, polinv, gamma;
    TEMPLATE(T, t) t;
    fmpz_t exp, q;
    int res = 0;
    slong i, k;

    if (pol->length <= 1)
    {
        TEMPLATE_PRINTF("Exception (%s_poly_factor_equal_deg_prob_vzgs): \n", T);
        flint_printf("Input polynomial is linear.\n");
        abort();
    }

    fmpz_init(q);
    TEMPLATE(T, ctx_order) (q, ctx);
    TEMPLATE(T, poly_init) (a, ctx);

    do
    {
        TEMPLATE(T, poly_randtest) (a, state, pol->length - 1, ctx);
    } while (a->length <= 1);

    TEMPLATE(T, poly_init) (polinv, ctx);
    TEMPLATE(T, poly_reverse) (polinv, pol, pol->length, ctx);
    TEMPLATE(T, poly_inv_series_newton) (polinv, polinv, polinv->length, ctx);

    
    TEMPLATE(T, poly_init) (b, ctx);
    TEMPLATE(T, poly_init) (xq, ctx);
    TEMPLATE(T, poly_init) (gamma, ctx);
    TEMPLATE(T, poly_powmod_xq_preinv)(xq, pol, polinv, ctx);
    TEMPLATE(T, poly_trace_frob_preinv)(b, a, d - 1, xq, pol, polinv, ctx);

    if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime) (ctx), 2) > 0)
    {
        /* compute a^{(q-1)/2} rem pol */
        fmpz_init(exp);
        fmpz_sub_ui(exp, q, 1);
        fmpz_fdiv_q_2exp(exp, exp, 1);

        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (gamma, b, exp, 0, pol,
                                                      polinv, ctx);

        /* Compute gcd(gamma, pol) */
        TEMPLATE(T, poly_gcd)(f, gamma, pol, ctx);
        if (f->length > 1 && f->length != pol->length)
        {
            res += 1;
            h2 = g;
        }
        else
            h2 = f;

        /* Compute gcd(gamma - 1, pol) */
        TEMPLATE(T, init) (t, ctx);
        TEMPLATE(T, sub_one) (t, gamma->coeffs, ctx);
        TEMPLATE(T, poly_set_coeff) (gamma, 0, t, ctx);
        TEMPLATE(T, clear) (t, ctx);

        TEMPLATE(T, poly_gcd)(h2, gamma, pol, ctx);
        if (h2->length > 1 && h2->length != pol->length)
            res += 1;
        
        fmpz_clear(exp);
    }
    else
    {
        /* compute gamma = (b^{2^{k-1}}+b^{2^{k-2}}+...+b^3+b^2+b) rem pol */
        k = TEMPLATE(T, ctx_degree) (ctx); 
        TEMPLATE(T, poly_init)(c, ctx);
        TEMPLATE(T, poly_set) (gamma, b, ctx);
        TEMPLATE(T, poly_set) (c, b, ctx);
        for (i = 1; i < k; i++)
        {
            /* c = b^{2^i} = (b^{2^{i-1}})^2 */
            TEMPLATE(T, poly_powmod_ui_binexp_preinv) (c, c, 2, pol, polinv,
                                                       ctx);
            TEMPLATE(T, poly_add) (gamma, b, c, ctx);
        }
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, poly_gcd)(f, gamma, pol, ctx);
        if (f->length > 1 && f->length != pol->length)
            res += 1;
    }

    TEMPLATE(T, poly_clear) (a, ctx);
    TEMPLATE(T, poly_clear) (b, ctx);
    TEMPLATE(T, poly_clear) (xq, ctx);
    TEMPLATE(T, poly_clear) (gamma, ctx);
    TEMPLATE(T, poly_clear) (polinv, ctx);
    fmpz_clear(q);

    return res;
}


#endif
