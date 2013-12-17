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

int
TEMPLATE(T, root)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op, slong n,
                  const TEMPLATE(T, ctx_t) ctx)
{
    int found = 0;
    slong i;
    fmpz_t q;
    TEMPLATE(T, poly_struct) * f;
    TEMPLATE(T, poly_t) v, finv, x, xqmx, g;
    TEMPLATE(T, poly_factor_t) sq_free, equal_deg;
    
    fmpz_init(q);
    TEMPLATE(T, ctx_order)(q, ctx);

    TEMPLATE(T, poly_init)(v, ctx);
    TEMPLATE(T, poly_init)(finv, ctx);
    TEMPLATE(T, poly_init)(x, ctx);
    TEMPLATE(T, poly_init)(xqmx, ctx);
    TEMPLATE(T, poly_init)(g, ctx);

    TEMPLATE(T, poly_fit_length)(v, n + 1, ctx);
    TEMPLATE(T, one)(v->coeffs + n, ctx);
    TEMPLATE(T, neg)(v->coeffs, op, ctx);
    _TEMPLATE(T, poly_set_length)(v, n + 1, ctx);
    
    /* Squarefree factorisation */
    TEMPLATE(T, poly_factor_init) (sq_free, ctx);
    TEMPLATE(T, poly_factor_squarefree) (sq_free, v, ctx);

    for (i = 0; i < sq_free->num; i++)
    {
        f = sq_free->poly + i;
        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_gen)(x, ctx);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(xqmx, x, q, 0,
                                                     f, finv, ctx);
        TEMPLATE(T, poly_sub)(xqmx, xqmx, x, ctx);

        TEMPLATE(T, poly_gcd)(g, xqmx, f, ctx);

        if (!TEMPLATE(T, poly_is_one)(g, ctx))
        {
            TEMPLATE(T, poly_factor_init)(equal_deg, ctx);
            TEMPLATE(T, poly_factor_equal_deg)(equal_deg, g, 1, ctx);
            if (equal_deg->num > 0)
            {
                TEMPLATE(T, neg)(rop, equal_deg->poly->coeffs, ctx);
                found = 1;
            }

            TEMPLATE(T, poly_factor_clear)(equal_deg, ctx);
            if (found)
                break;
        }
    }

    TEMPLATE(T, poly_factor_clear)(sq_free, ctx);
    TEMPLATE(T, poly_clear)(x, ctx);
    TEMPLATE(T, poly_clear)(xqmx, ctx);
    TEMPLATE(T, poly_clear)(v, ctx);
    TEMPLATE(T, poly_clear)(finv, ctx);
    TEMPLATE(T, poly_clear)(g, ctx);
    fmpz_clear(q);
    
    return found;
}


#endif
